#include <torch/all.h>
#include <torch/script.h>
#include <cassert>
#include <torch/nn/module.h>
#include <torch/nn/functional.h>
#include <torch/serialize.h>
#include "c_libtorch.h"
#include "ANIModel.h"
//#include <torch/torch.h>
#include <torch/serialize/archive.h>
#include <torch/serialize/tensor.h>
#include <iostream>
#include <vector>
#include <memory>
#include <chrono>
#include <cstdint>
#include <omp.h>

using namespace torch::autograd;


struct Atom {
    std::string species;
    std::vector<float> coordinates;
};


void ANIModel::load_model(int model_type, const std::string &nn_path) {
    if (model_type == 1) {
        module = torch::jit::load(nn_path + "ANI1x_saved2.pt");
	std::cout << "Model loaded from: " << nn_path + "ANI1x_saved2.pt" << std::endl;
    } else if (model_type == 2) {
        module = torch::jit::load(nn_path + "ANI2x_saved.pt");
	std::cout << "Model loaded from: " << nn_path + "ANI2x_saved.pt" << std::endl;
    } else {
        std::cerr << "Invalid model type!" << std::endl;
    }
}

void ANIModel::load_custom_model(const std::string &aev_name, const std::string &model_name, const std::string &nn_path) {
        aev_computer = torch::jit::load(nn_path + aev_name);
        module = torch::jit::load(nn_path + model_name);
	std::cout << "AEV loaded from: " << nn_path + aev_name << std::endl;
	std::cout << "Custom model loaded from: " << nn_path + model_name << std::endl;
}

/*
 // using back propagation
void ANIModel::get_energy_grad(const torch::Tensor& coordinates, 
                               const torch::Tensor& species, 
                               float* atomic_energies, 
                               float* gradients, 
                               float* forces, 
                               int num_atoms) {
    std::vector<torch::jit::IValue> inputs;
    inputs.push_back(std::make_tuple(species, coordinates));

    auto output = module.forward(inputs).toTuple();
    at::Tensor energy_tensor = output->elements()[1].toTensor();

    auto energy = energy_tensor.item<float>();
    energy_tensor.backward(torch::ones_like(energy_tensor));

    auto gradient = coordinates.grad();

    if (!gradient.defined() || gradient.numel() == 0) {
        std::cerr << "Error: Gradient is not defined or empty." << std::endl;
        return;
    }

    auto force = -gradient;
    auto atomic_energies_tensor = module.get_method("atomic_energies")(inputs).toTuple()->elements()[1].toTensor();

    std::cout << "=========TESTING FOR OBJECT BASED MODEL LOADING ===============" << std::endl;
    std::cout << " Energy: " << energy << std::endl;
    std::cout << " Force: " << force << std::endl;

    memcpy(atomic_energies, atomic_energies_tensor.data_ptr<float>(), atomic_energies_tensor.numel() * sizeof(float));
    memcpy(gradients, gradient.data_ptr<float>(), gradient.numel() * sizeof(float));
    memcpy(forces, force.data_ptr<float>(), force.numel() * sizeof(float));
    coordinates.grad().zero_();
}
*/

// using autograd
void ANIModel::get_energy_grad(const torch::Tensor& coordinates,
                               const torch::Tensor& species,
                               float* atomic_energies,
                               float* gradients,
                               float* forces,
                               int num_atoms) {

//  auto start_total = std::chrono::high_resolution_clock::now(); // Start total timer
//  Forward pass
//  auto start_forward = std::chrono::high_resolution_clock::now();
    std::vector<torch::jit::IValue> inputs;
    inputs.push_back(std::make_tuple(species, coordinates));

    auto output = module.forward(inputs).toTuple();
    at::Tensor energy_tensor = output->elements()[1].toTensor();
//  auto end_forward = std::chrono::high_resolution_clock::now();
//  std::chrono::duration<double> forward_duration = end_forward - start_forward;
//  std::cout << "Forward pass time: " << forward_duration.count() << " seconds" << std::endl;

//  integrate GPU for the autograd part
//  auto start_grad = std::chrono::high_resolution_clock::now();

    auto  gradients_vec = torch::autograd::grad({energy_tensor}, {coordinates}, {}, true, false);
    torch::Tensor gradient = gradients_vec[0];

//  auto end_grad = std::chrono::high_resolution_clock::now();
//  std::chrono::duration<double> grad_duration = end_grad - start_grad;
//  std::cout << "Gradient calculation time: " << grad_duration.count() << " seconds" << std::endl;

    if (!gradient.defined() || gradient.numel() == 0) {
        std::cerr << "Error: Gradient is not defined or empty." << std::endl;
        return;
    }

//  Force calcn
//  auto start_force = std::chrono::high_resolution_clock::now();
    torch::Tensor force = -gradient;
//  auto end_force = std::chrono::high_resolution_clock::now();
//  std::chrono::duration<double> force_duration = end_force - start_force;
//  std::cout << "Force calculation time: " << force_duration.count() << " seconds" << std::endl;

//  Measure the time for atomic energy retrieval
//  auto start_atomic_energy = std::chrono::high_resolution_clock::now();        
    
    auto atomic_energies_tensor = module.get_method("atomic_energies")(inputs).toTuple()->elements()[1].toTensor();

//  auto end_atomic_energy = std::chrono::high_resolution_clock::now();
//  std::chrono::duration<double> atomic_energy_duration = end_atomic_energy - start_atomic_energy;
//  std::cout << "Atomic energy retrieval time: " << atomic_energy_duration.count() << " seconds" << std::endl;

//  Print total time
//  auto end_total = std::chrono::high_resolution_clock::now();
//  std::chrono::duration<double> total_duration = end_total - start_total;
//  std::cout << "Total time for get_energy_grad: " << total_duration.count() << " seconds" << std::endl;

    std::cout << "=========TESTING FOR OBJECT BASED MODEL LOADING ===============" << std::endl;
    std::cout << " Energy: " << energy_tensor.item<float>() << std::endl;
    std::cout << " Force: " << force << std::endl;
    std::cout << "=========================================" << std::endl;

    memcpy(atomic_energies, atomic_energies_tensor.data_ptr<float>(), atomic_energies_tensor.numel() * sizeof(float));
    memcpy(gradients, gradient.data_ptr<float>(), gradient.numel() * sizeof(float));
    memcpy(forces, force.data_ptr<float>(), force.numel() * sizeof(float));
}

void ANIModel::get_custom_energy_grad(float* coordinates_data, int64_t* species_data, float* elecpots_data, int num_atoms, float* custom_energy, float* cus_grads, float* cus_forces) {

    torch::Tensor coordinates = torch::from_blob(coordinates_data, {1, num_atoms, 3}, torch::kFloat32).clone().set_requires_grad(true);
    torch::Tensor species = torch::from_blob(species_data, {1, num_atoms}, torch::kInt64).clone();
    torch::Tensor elecpots = torch::from_blob(elecpots_data, {1, num_atoms}, torch::kFloat32).clone();

    coordinates = coordinates.contiguous();

    std::map<int, double> ani1x_sae_dict_byIdx = {
                {0, -0.60095298}, // H
                {1, -38.08316124}, // C
                {2, -54.7077577}, // N
                {3, -75.19446356} // O
        };

    double shift = 0.0;
    for (int i = 0; i < species.size(1); ++i) {
        int atom_type = species[0][i].item<int>();
        shift += ani1x_sae_dict_byIdx[atom_type];
    }	

    auto aev_input = std::make_tuple(species, coordinates);
    auto aev_output = aev_computer.forward({aev_input}).toTuple();
    torch::Tensor aevs = aev_output->elements()[1].toTensor();  // Get AEV output

    torch::Tensor aep = torch::cat({aevs, elecpots.unsqueeze(-1)}, -1);

    auto model_input = std::make_tuple(species, aep);
    auto energy_output = module.forward({model_input}).toTuple();
 
    torch::Tensor energy_unshifted = energy_output->elements()[1].toTensor(); //c
    torch::Tensor energy_shifted = energy_unshifted + shift; //c

    std::cout << "Energy (unshifted): " << energy_unshifted.item<float>() << std::endl; //c
    std::cout << "Energy (shifted): " << energy_shifted.item<float>() << std::endl; //c

    std::vector<torch::Tensor> gradients = torch::autograd::grad({energy_shifted}, {coordinates});
    torch::Tensor derivative = gradients[0];

    torch::Tensor force = -derivative;

    std::cout << "Force: " << force << std::endl; //c

    memcpy(cus_grads, derivative.data_ptr<float>(), derivative.numel() * sizeof(float));
    memcpy(cus_forces, force.data_ptr<float>(), force.numel() * sizeof(float));
    memcpy(custom_energy, energy_shifted.data_ptr<float>(), sizeof(float));

}

// Hardcoded custom model routine 
void engrad_custom_model(float* coordinates_data, int64_t* species_data, float* elecpots_data, int num_atoms, float* custom_energy, float* cus_grads, float* cus_forces) {

    torch::jit::script::Module aev_computer = torch::jit::load("/depot/lslipche/data/skp/torch_skp_branch/libefp/nnlib/aev_scripted.pt");
    torch::jit::script::Module model = torch::jit::load("/depot/lslipche/data/skp/torch_skp_branch/libefp/nnlib/custom_model_script.pt");

    torch::Tensor coordinates = torch::from_blob(coordinates_data, {1, num_atoms, 3}, torch::kFloat32).clone().set_requires_grad(true);
    torch::Tensor species = torch::from_blob(species_data, {1, num_atoms}, torch::kInt64).clone();
    torch::Tensor elecpots = torch::from_blob(elecpots_data, {1, num_atoms}, torch::kFloat32).clone();

    coordinates = coordinates.contiguous();

    auto aev_input = std::make_tuple(species, coordinates);
    auto aev_output = aev_computer.forward({aev_input}).toTuple();
    torch::Tensor aevs = aev_output->elements()[1].toTensor();  // Get AEV output

    torch::Tensor aep = torch::cat({aevs, elecpots.unsqueeze(-1)}, -1);

    auto model_input = std::make_tuple(species, aep);
    auto energy_output = model.forward({model_input}).toTuple();
    torch::Tensor energy = energy_output->elements()[1].toTensor();

    std::vector<torch::Tensor> gradients = torch::autograd::grad({energy}, {coordinates});
    torch::Tensor derivative = gradients[0];

    torch::Tensor force = -derivative;

    memcpy(cus_grads, derivative.data_ptr<float>(), derivative.numel() * sizeof(float));
    memcpy(cus_forces, force.data_ptr<float>(), force.numel() * sizeof(float));
    memcpy(custom_energy, energy.data_ptr<float>(), sizeof(float));

}


torch::jit::script::Module loadModel(const std::string& modelPath) {
    try {

	// Load the model from file
	torch::jit::script::Module model;
        model = torch::jit::load(modelPath);
        model.to(torch::kCPU);
        model.eval();
        return model;
    } catch (const c10::Error& e) {
        std::cerr << "Error loading the model: " << e.what() << std::endl;
        exit(1);
    }
}

void generateEnergyForces(torch::jit::script::Module& model, const std::vector<std::vector<float>>& coordinates,
                          float& energy, std::vector<std::vector<float>>& forces) {
    try {
        std::vector<std::vector<float>> coordCopy(coordinates);
        torch::Tensor inputTensor = torch::from_blob(coordCopy.data(), {1, static_cast<long>(coordCopy.size()), 3});
	std::vector<c10::IValue> energyInputs = {inputTensor};
        torch::Tensor energyTensor = model.forward(energyInputs).toTensor();

        // Calculate forces using autograd
	std::vector<torch::Tensor> forceInputs = {energyTensor};
        std::vector<torch::Tensor> forcesTensor = torch::autograd::grad(forceInputs, {inputTensor}, {torch::ones_like(energyTensor)}, true, true);

        // Extract energy and forces from tensors
	energy = energyTensor.item<float>();
        forces.resize(coordinates.size());
        for (size_t i = 0; i < coordinates.size(); ++i) {
            forces[i].resize(3);
            for (size_t j = 0; j < 3; ++j) {
                forces[i][j] = forcesTensor[0][i][j].item<float>();
            }
        }
    } catch (const c10::Error& e) {
        std::cerr << "Error generating energy and forces: " << e.what() << std::endl;
        exit(1);
    }
}

void generateSpeciesEnergyForces(torch::jit::script::Module& model,
                          const torch::Tensor& species,
                          const torch::Tensor& coordinates,
                          float& energy,
                          torch::Tensor& forces) {
    try {
        torch::IValue species_coordinates = std::make_tuple(species, coordinates);
	auto result = model.forward({species_coordinates}).toTuple();
	energy = result->elements()[0].toTensor().index({0}).item<float>();  // Extract the scalar value
	forces = result->elements()[1].toTensor().clone();  // Make a clone to avoid modifying the original tensor
	
  
    } catch (const c10::Error& e) {
        std::cerr << "Error generating energy and forces: " << e.what() << std::endl;
        exit(1);
    }
}


int64_t mapSpeciesToInteger(const std::string& species) {
    if (species == "H")
        return 1;
    else if (species == "O")
        return 2;
    else if (species == "C")
        return 3;
    else
        return 0; 
}

struct Net : torch::nn::Module {
    Net() {
        fc1 = register_module("fc1", torch::nn::Linear(3, 16));
        fc2 = register_module("fc2", torch::nn::Linear(16, 3));
    }

    torch::Tensor forward(torch::Tensor x) {
        x = torch::relu(fc1->forward(x));
        x = fc2->forward(x);
        return x;
    }

    torch::nn::Linear fc1{nullptr}, fc2{nullptr};
};

void trainModel(Net &model, const std::vector<std::vector<float>> &input_data, const std::vector<float> &target_data, int num_epochs, float learning_rate) {
    torch::optim::SGD optimizer(model.parameters(), torch::optim::SGDOptions(learning_rate));

    for (int epoch = 0; epoch < num_epochs; ++epoch) {
        std::vector<std::vector<float>> non_const_input_data = input_data;
        torch::Tensor inputs = torch::from_blob(non_const_input_data.data(), {static_cast<long>(non_const_input_data.size()), static_cast<long>(non_const_input_data[0].size())});

        std::vector<float> target_data_nonconst(target_data.begin(), target_data.end());
        torch::TensorOptions options = torch::TensorOptions().dtype(torch::kFloat32);
        torch::Tensor targets = torch::from_blob(target_data_nonconst.data(), {static_cast<long>(target_data_nonconst.size())}, options);
        optimizer.zero_grad();
        torch::Tensor output = model.forward(inputs);
        torch::Tensor loss = torch::mse_loss(output, targets);
        loss.backward();
        optimizer.step();
        if (epoch % 10 == 0) {
            std::cout << "Epoch: " << epoch << ", Loss: " << loss.item<float>() << std::endl;
        }
    }
}


void generateEnergy(Net &model, const std::vector<std::vector<float>> &input_data) {
    torch::NoGradGuard no_grad;

    int batch_size = input_data.size(); // del later
    int input_size = input_data[0].size(); // del later

    std::vector<std::vector<float>> non_const_input_data = input_data;

    std::vector<float> flattened_input_data;
    flattened_input_data.reserve(batch_size * input_size);
    for (const auto &sample : non_const_input_data) {
        flattened_input_data.insert(flattened_input_data.end(), sample.begin(), sample.end());
    }


   torch::Tensor inputs = torch::from_blob(flattened_input_data.data(), {batch_size, input_size}); // de later



    torch::Tensor output = model.forward(inputs);

    std::cout << "Energy: ";
    std::vector<float> energy_values(output.data_ptr<float>(), output.data_ptr<float>() + output.numel());

    for (const auto &value : energy_values) {
        std::cout << value << " ";
    }
    std::cout << std::endl;
}


torch::Tensor compute_gradient(const torch::Tensor& input) {
     // Enable gradient computation
     torch::autograd::GradMode::set_enabled(true);
     // Create a variable from the input tensor
     torch::Tensor variable = input.clone().detach().requires_grad_(true);
     // Compute the output tensor
     torch::Tensor output = variable * variable/* your computation using the variable */;
     // Compute gradients of the output tensor with respect to the variable
     torch::autograd::variable_list grad_outputs = {torch::ones_like(output)};
     torch::autograd::variable_list gradients = torch::autograd::grad({output}, {variable}, grad_outputs, /* retain_graph */ true);
     torch::Tensor gradient = gradients[0];
     
     return gradient;
}     


// Hardcoded ANI1 routine                                            
void get_ANI1_energy_grad(const torch::Tensor& coordinates, const torch::Tensor& species, float *atomic_energies, float *gradients, float *forces) {
  
	//const char* model_path_env = std::getenv("TORCHANI_DIR");
        //std::string model_path = std::string(model_path_env) + "ANI1x_saved2.pt";
        //std::cout << "Model loaded successfully from " << model_path << std::endl;

//    torch::jit::Module module = torch::jit::load(model_path);
    torch::jit::Module module = torch::jit::load("/depot/lslipche/data/skp/torch_skp_branch/libefp/nnlib/ANI1x_saved2.pt");
    std::vector<torch::jit::IValue> inputs;
    inputs.push_back(std::make_tuple(species, coordinates));

    auto output = module.forward(inputs).toTuple();
    at::Tensor energy_tensor = output->elements()[1].toTensor();

    auto energy = energy_tensor.item<float>();
    energy_tensor.backward(torch::ones_like(energy_tensor));

    auto gradient = coordinates.grad();

    if (!gradient.defined() || gradient.numel() == 0) {
        std::cerr << "Error: Gradient is not defined or empty." << std::endl;
        return;
    }

    auto force = -gradient;
    auto atomic_energies_tensor = module.get_method("atomic_energies")(inputs).toTuple()->elements()[1].toTensor();

    memcpy(atomic_energies, atomic_energies_tensor.data_ptr<float>(), atomic_energies_tensor.numel() * sizeof(float));
    memcpy(gradients, gradient.data_ptr<float>(), gradient.numel() * sizeof(float));
    memcpy(forces, force.data_ptr<float>(), force.numel() * sizeof(float));
    coordinates.grad().zero_();

}

// Hardcoded ANI2 routine
void get_ANI2_energy_grad(const torch::Tensor& coordinates, const torch::Tensor& species, float *atomic_energies, float *gradients, float *forces) {

    torch::jit::Module module = torch::jit::load("/depot/lslipche/data/skp/torch_skp_branch/libefp/nnlib/ANI2x_saved.pt");
    std::vector<torch::jit::IValue> inputs;
    inputs.push_back(std::make_tuple(species, coordinates));

    auto output = module.forward(inputs).toTuple();
    at::Tensor energy_tensor = output->elements()[1].toTensor();

    auto energy = energy_tensor.item<float>();
    energy_tensor.backward(torch::ones_like(energy_tensor));

    auto gradient = coordinates.grad();

    if (!gradient.defined() || gradient.numel() == 0) {
        std::cerr << "Error: Gradient is not defined or empty." << std::endl;
        return;
    }

    auto force = -gradient;
    auto atomic_energies_tensor = module.get_method("atomic_energies")(inputs).toTuple()->elements()[1].toTensor();

    memcpy(atomic_energies, atomic_energies_tensor.data_ptr<float>(), atomic_energies_tensor.numel() * sizeof(float));
    memcpy(gradients, gradient.data_ptr<float>(), gradient.numel() * sizeof(float));
    memcpy(forces, force.data_ptr<float>(), force.numel() * sizeof(float));
    coordinates.grad().zero_();
 
}


//================================= END OF C++ routines ==================================//
 

extern "C" {

struct TensorData
{
  at::Tensor tensor;
};

// SKP wrappers
 
void *compute_gradient_c(float* data, int64_t* sizes, int ndim) {

  torch::TensorOptions options;
  options = options.dtype(torch::kFloat32);
  torch::Tensor input_tensor = torch::from_blob(data, torch::IntArrayRef(sizes, ndim), options).clone();

  torch::Tensor gradient_tensor = compute_gradient(input_tensor);

  size_t tensor_size = gradient_tensor.numel();
        float* gradient_data = new float[tensor_size];
        std::memcpy(gradient_data, gradient_tensor.data_ptr<float>(), tensor_size * sizeof(float));

  int64_t* gradient_sizes = new int64_t[gradient_tensor.dim()];
        std::memcpy(gradient_sizes, gradient_tensor.sizes().data(), gradient_tensor.dim() * sizeof(int64_t));

        struct Tensor {
            void* data;
            int64_t* sizes;
            int ndim;
            int type_id;
            int is_variable;
        }; 

   Tensor* gradient = new Tensor();
        gradient->data = gradient_data;
        gradient->sizes = gradient_sizes;
        gradient->ndim = gradient_tensor.dim();
        gradient->type_id = static_cast<int>(gradient_tensor.scalar_type());
        gradient->is_variable = 0;

        return gradient;
}

void destroy_tensor(struct Tensor *tensor) {
        //delete[] tensor->data;
	delete[] static_cast<float*>(tensor->data);
        delete[] tensor->sizes;
        delete tensor;
}

//  SKP started on June 5th=========================//

Net *createNet() {
    return new Net();
}

void destroyNet(Net *model) {
    delete model;
}

void forward(Net *model, const float *inputs, float *output, int input_size, int output_size) {
    torch::Tensor input_tensor = torch::from_blob(const_cast<float*>(inputs), {input_size});
    torch::Tensor output_tensor = torch::from_blob(output, {output_size});
    torch::Tensor result = model->forward(input_tensor);

    std::memcpy(output, result.data_ptr<float>(), output_size * sizeof(float));
}

void trainModelWrapper(Net *model, const float **input_data, const float *target_data,int num_samples, int num_epochs, float learning_rate) {

    std::vector<std::vector<float>> input_vec;
    for (int i = 0; i < num_samples; ++i) {
        std::vector<float> data;
        for (int j = 0; j < 3; ++j) {
            data.push_back(input_data[i][j]);
        }
        input_vec.push_back(data);
    }

    trainModel(*model, input_vec, std::vector<float>(target_data, target_data + num_samples), num_epochs, learning_rate);

}

void generateEnergyWrapper(Net* model, const float **input_data, int batch_size, int input_size) {
    std::vector<std::vector<float>> input_vec(batch_size, std::vector<float>(input_size));

    for (int i = 0; i < batch_size; ++i) {
            for (int j = 0; j < input_size; ++j) {
                input_vec[i][j] = input_data[i][j];
            }
        }

    generateEnergy(*model, input_vec);
}

// June6th

void *loadModelWrapper(const char *modelPath) {
    std::string path(modelPath);
    torch::jit::script::Module *model = new torch::jit::script::Module(loadModel(path));
    return static_cast<void*>(model);
}

void generateEnergyForcesWrapper(const void* model, const float* const* coordinates, int num_atoms, float* energy, float* const* forces) {

    const torch::jit::script::Module* torchModel = static_cast<const torch::jit::script::Module*>(model);

    std::vector<std::vector<float>> coordVec(num_atoms, std::vector<float>(3));
    for (int i = 0; i < num_atoms; ++i) {
        for (int j = 0; j < 3; ++j) {
            coordVec[i][j] = coordinates[i][j];
        }
    }

    std::vector<std::vector<float>> forcesVec;
    generateEnergyForces(const_cast<torch::jit::script::Module&>(*torchModel), coordVec, *energy, forcesVec);

    for (int i = 0; i < num_atoms; ++i) {
        for (int j = 0; j < 3; ++j) {
            forces[i][j] = forcesVec[i][j];
        }
    }
}

void generateSpeciesEnergyForcesWrapper(const void* model,
                                 const float* const* coordinates,
                                 const int* species,
                                 int num_atoms,
                                 float* energy,
                                 float* const* forces) {
    const torch::jit::script::Module* torchModel = static_cast<const torch::jit::script::Module*>(model);

    // Convert the species and coordinates to Torch tensors
//    torch::Tensor speciesTensor = torch::from_blob(species, {1, num_atoms}, torch::kInt32).clone();
//    torch::Tensor coordinatesTensor = torch::from_blob(coordinates, {1, num_atoms, 3}).clone();

    torch::Tensor speciesTensor = torch::from_blob(const_cast<int*>(species), {1, num_atoms}, torch::kInt32).clone();
    torch::Tensor coordinatesTensor = torch::from_blob(const_cast<float**>(coordinates), {1, num_atoms, 3}).clone();

    torch::jit::script::Module& nonConstModel = const_cast<torch::jit::script::Module&>(*torchModel);

    // Generate energy and forces
    torch::Tensor forcesTensor;
    generateSpeciesEnergyForces(nonConstModel, speciesTensor, coordinatesTensor, *energy, forcesTensor);

    // Copy the forces to the output array
    for (int i = 0; i < num_atoms; ++i) {
        for (int j = 0; j < 3; ++j) {
//            forces[i][j] = forcesTensor[0][i][j].item<float>();
	   forces[i][j] = forcesTensor[i][j].item<float>();
        }
    }
}

//======== SKP June 29 =================================//

// previously nnp_test7_wrapper
// this routine should have module in argument rather than model type
// and probably get_ANI1_energy_grad and get_ANI2_energy_grad will boil down to just get_ANI_energy_grad(.., .., ..., module)
// load NNP and its wrapper.. call that wrapper in opt.c/main.c...

void get_torch_energy_grad(float* coordinates_data, int* species_data, int num_atoms, float *atomic_energies, float *gradients, float *forces, int model_type) {

	torch::Tensor speciesTensor = torch::from_blob(const_cast<int*>(species_data), {1, num_atoms}, torch::kInt32);
        torch::Tensor coordinatesTensor = torch::from_blob(const_cast<float*>(coordinates_data), {1, num_atoms, 3}, torch::requires_grad(true));

	if(model_type == 1) get_ANI1_energy_grad(coordinatesTensor, speciesTensor, atomic_energies, gradients, forces);
	if(model_type == 2) get_ANI2_energy_grad(coordinatesTensor, speciesTensor, atomic_energies, gradients, forces);

}


void engrad_custom_model_wrapper(float* coordinates_data, int64_t* species_data, float* elecpots_data, int num_atoms, float* custom_energy, float* gradients, float* forces) {

    engrad_custom_model(coordinates_data, species_data, elecpots_data, num_atoms, custom_energy, gradients, forces);
//    std::cout << "Custom energy in wrapper " << *custom_energy << std::endl;
}

//=================================================


ANIModel* ANIModel_new() {
    return new ANIModel();
}


void load_ani_model(ANIModel* model, int model_type, const char* nn_path) {
    std::string path_str(nn_path);
    model->load_model(model_type, path_str);
}

void load_custom_ani_model(ANIModel* model, const char* aev_name, const char* model_name, const char* nn_path) {
    std::string path_str(nn_path);
    std::string mod_str(model_name);
    std::string aev_str(aev_name);

    model->load_custom_model(aev_str, mod_str, path_str);
}


void get_ani_energy_grad(ANIModel* model, float* coordinates, int* species, float* atomic_energies, float* gradients, float* forces, int num_atoms) {
    auto coordinates_tensor = torch::from_blob((float*)coordinates, {1, num_atoms, 3}, torch::requires_grad(true));
    auto species_tensor = torch::from_blob((int*)species, {1, num_atoms}, torch::kInt32);

    auto start_ani = std::chrono::high_resolution_clock::now();
 
    model->get_energy_grad(coordinates_tensor, species_tensor, atomic_energies, gradients, forces, num_atoms);

    auto end_ani = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> ani_duration = end_ani - start_ani;
    std::cout << "get_energy_grad calculation time: " << ani_duration.count() << " seconds" << std::endl;

    std::cout << "==========END OF TEST OBJECT BASED MODEL LOADING=============" << std::endl; 
}
 
//void get_ani_energy_grad(ANIModel* model, float* coordinates, int64_t* species, float* ani_energy, float* gradients, float* forces, int num_atoms) {
 
//    model->get_energy_grad(coordinates, species, ani_energy, gradients, forces, num_atoms);

//}
 
void get_custom_energy_grad_wrapper(ANIModel* model, float* coordinates, int64_t* species, float* elecpots, int num_atoms, float* custom_energy, float* gradients, float* forces) {

    auto start_custom = std::chrono::high_resolution_clock::now();

    model->get_custom_energy_grad(coordinates, species, elecpots, num_atoms, custom_energy, gradients, forces); 
 
    auto end_custom = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> custom_duration = end_custom - start_custom;
    std::cout << "get_custom_energy_grad calculation time: " << custom_duration.count() << " seconds" << std::endl;    
}

void ANIModel_delete(ANIModel* model) {
    delete model;
}


//=====================================// 



} // extern "C"
