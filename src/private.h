/*-
 * Copyright (c) 2012-2013 Ilya Kaliman
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 */

#ifndef LIBEFP_PRIVATE_H
#define LIBEFP_PRIVATE_H

#include <assert.h>

#include "bvec.h"
#include "efp.h"
#include "ff.h"
#include "int.h"
#include "log.h"
#include "swf.h"
#include "terms.h"
#include "util.h"

#define EFP_EXPORT __attribute__((visibility("default")))
#define ARRAY_SIZE(arr) (sizeof(arr)/sizeof(arr[0]))

struct frag {
	/* fragment name */
	char name[32];

	/* fragment center of mass */
	double x, y, z;

	/* rotation matrix representing orientation of a fragment */
	mat_t rotmat;

	/* pointer to the initial fragment state in library */
	const struct frag *lib;

	/* force on fragment center of mass */
	vec_t force;

	/* torque on fragment */
	vec_t torque;

	/* number of atoms in this fragment */
	size_t n_atoms;

	/* fragment atoms */
	struct efp_atom *atoms;

	/* distributed multipoles */
	struct multipole_pt {
		double x, y, z;
		double monopole;
		vec_t dipole;
		double quadrupole[6];
		double octupole[10];
	} *multipole_pts;

	/* number of distributed multipole points */
	size_t n_multipole_pts;

	/* electrostatic screening parameters */
	double *screen_params;

	/* ab initio electrostatic screening parameters */
	double *ai_screen_params;

	/* distributed polarizability points */
	struct polarizable_pt {
		double x, y, z;
		mat_t tensor;
		vec_t elec_field;
		vec_t elec_field_wf;
		vec_t induced_dipole;
		vec_t induced_dipole_conj;
	} *polarizable_pts;

	/* number of distributed polarizability points */
	size_t n_polarizable_pts;

	/* dynamic polarizability points */
	struct dynamic_polarizable_pt {
		double x, y, z;
		double trace[12];
	} *dynamic_polarizable_pts;

	/* number of dynamic polarizability points */
	size_t n_dynamic_polarizable_pts;

	/* number of localized molecular orbitals */
	size_t n_lmo;

	/* localized molecular orbital centroids */
	vec_t *lmo_centroids;

	/* spin multiplicity */
	int multiplicity;

	/* number of exchange repulsion basis shells */
	size_t n_xr_shells;

	/* exchange repulsion basis shells */
	struct shell *xr_shells;

	/* upper triangle of fock matrix, size = n_lmo * (n_lmo + 1) / 2 */
	double *xr_fock_mat;

	/* exchange repulsion wavefunction size */
	size_t xr_wf_size;

	/* exchange repulsion wavefunction, size = n_lmo * xr_wf_size */
	double *xr_wf;

	/* rotational derivatives of MO coefficients */
	double *xr_wf_deriv[3];

	/* offset of polarizable points for this fragment */
	size_t polarizable_offset;

	struct ff_atom {
		char type[32]; /* atom type in force field */
		size_t idx;    /* index in atoms array */
	} *ff_atoms;

	size_t n_ff_atoms;

	struct ff_link {
		size_t idx1;   /* index in ff_atoms array */
		size_t idx2;   /* index in ff_atoms array */
	} *ff_links;

	size_t n_ff_links;

	/* offset in array of ff atoms */
	size_t ff_offset;
};

struct efp {
	/* number of fragments */
	size_t n_frag;

	/* array of fragments */
	struct frag *frags;

	/* number of fragments in the library */
	size_t n_lib;

	/* array with the library of fragment initial parameters */
	struct frag **lib;

	/* callback which computes electric field from electrons */
	efp_electron_density_field_fn get_electron_density_field;

	/* user data for get_electron_density_field */
	void *get_electron_density_field_user_data;

	/* user parameters for this EFP computation */
	struct efp_opts opts;

	/* gradient will also be computed if nonzero */
	int do_gradient;

	/* periodic simulation box size */
	vec_t box;

	/* stress tensor */
	mat_t stress;

	/* number of point charges */
	size_t n_ptc;

	struct point_charge {
		double x, y, z;
		double charge;
		vec_t grad;
	} *point_charges;

	/* total number of polarizable points */
	size_t n_polarizable_pts;

	/* EFP energy terms */
	struct efp_energy energy;

	/* force field for covalent links */
	struct ff *ff;

	/* bit matrix of covalent links between fragments */
	struct bvec *links_bvec;
};

#endif /* LIBEFP_PRIVATE_H */
