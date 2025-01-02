/*-
 * Copyright (c) 2012-2017 Ilya Kaliman
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
#include <stdbool.h>

#include "efp.h"
#include "int.h"
#include "log.h"
#include "swf.h"
#include "terms.h"
#include "util.h"

#define EFP_EXPORT

#define ARRAY_SIZE(arr) (sizeof(arr)/sizeof(arr[0]))

struct multipole_pt {
    char label[32];   /**< Multipole point label. */
    double x, y, z;
    double znuc;
	double monopole;
	vec_t dipole;
	double quadrupole[6];
	double octupole[10];
	double screen2;   /**<Electrostatic frag-frag screening parameter>*/
	double screen0;    /**<Electrostatic AI-frag screening parameter>*/
    bool if_znuc;
    bool if_mon;
    bool if_dip;
    bool if_quad;
    bool if_oct;
    bool if_scr2;
    bool if_scr0;
};

struct polarizable_pt {
	double x, y, z;
	mat_t tensor;
	vec_t elec_field;
	vec_t elec_field_wf;  // electric field due to wavefunction of the QM region in the current state
	vec_t ligand_field;   // field due to ligand
	vec_t indip;          // current induced dipole for the current electronic state
	vec_t indip_old;      // previous induced dipole for the current electronic state
	vec_t indipconj;
	vec_t indipconj_old;
    vec_t indip_gs;      // induced dipole for the ground electronic state
    vec_t indipconj_gs;  // induced dipole for the ground electronic state
};

/* polarizable point on a ligand to store fields due to other fragments */
struct ligand_pt {
    /* fields of all (n_frag) fragments on this point  */
    vec_t *fragment_field;

    /* number of fragments producing fields */
    size_t n_frag;
};

struct dynamic_polarizable_pt {
	double x, y, z;
	mat_t tensor[12];
};

struct dipquad_polarizable_pt {
    double x, y, z;
    t3_t tensor[12];
};

struct ff_atom {
	char type[32]; /* atom type in force field */
	size_t idx;    /* index in atoms array */
};

struct ff_link {
	size_t idx1;   /* index in ff_atoms array */
	size_t idx2;   /* index in ff_atoms array */
};

struct frag {
	/* fragment name */
	char name[32];

	/* fragment center of mass */
	double x, y, z;

	/* rotation matrix representing orientation of a fragment */
	mat_t rotmat;

	/* pointer to the initial fragment state in library */
	const struct frag *lib;

    /* pointer to the updated fragment parameters */
    const struct frag *lib_current;

    /* rmsd between initial and updated structures of the fragment */
    double rmsd;

    /* number of atoms in this fragment */
	size_t n_atoms;

	/* fragment atoms */
	struct efp_atom *atoms;

	/* distributed multipoles */
	struct multipole_pt *multipole_pts;

	/* number of distributed multipole points */
	size_t n_multipole_pts;

	/* polarization damping parameter */
	double pol_damp;

	/* distributed polarizability points */
	struct polarizable_pt *polarizable_pts;

	/* number of distributed polarizability points */
	size_t n_polarizable_pts;

    /* dynamic polarizability points */
	struct dynamic_polarizable_pt *dynamic_polarizable_pts;

    /* dipole-quadrupole dynamic polarizability points */
    struct dipquad_polarizable_pt *dipquad_polarizable_pts;

    /* number of dynamic polarizability points */
	size_t n_dynamic_polarizable_pts;

    /* number of dipole-quadrupole dynamic polarizability points */
    size_t n_dipquad_polarizable_pts;

    /* number of localized molecular orbitals */
	size_t n_lmo;

	/* localized molecular orbital centroids */
	vec_t *lmo_centroids;

	/* spin multiplicity */
	int multiplicity;

	/* number of exchange repulsion atoms */
	size_t n_xr_atoms;

	/* exchange repulsion atoms */
	struct xr_atom *xr_atoms;

	/* upper triangle of fock matrix, size = n_lmo * (n_lmo + 1) / 2 */
	double *xr_fock_mat;

	/* exchange repulsion wavefunction size */
	size_t xr_wf_size;

	/* exchange repulsion wavefunction, size = n_lmo * xr_wf_size */
	double *xr_wf;

	/* rotational derivatives of MO coefficients */
	double *xr_wf_deriv[3];

	/* fitted ai-efp exchange-repulsion parameters */
	double *xrfit;

	/* offset of polarizable points for this fragment */
	size_t polarizable_offset;

};

/* structure derived from struct frag for describing ligand */
struct ligand {
    /* fragment */
    struct frag *ligand_frag;

    /* array of ligand points */
    struct ligand_pt *ligand_pts;

    /* number of ligand points */
    size_t n_ligand_pts;
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

    /* number of current (updated) fragments in the library */
    size_t n_lib_current;

    /* array with the library of fragment updated (shifted) parameters */
    struct frag **lib_current;

    /* pointer to ligand fragment */
    struct ligand *ligand;

    /* ligand index in fragment list */
    size_t ligand_index;
    
    /* callback which computes electric field from electrons */
    efp_electron_density_field_fn get_electron_density_field;

    /* user data for get_electron_density_field */
    void *get_electron_density_field_user_data;

    /* user parameters for this EFP computation */
	struct efp_opts opts;

	/* gradient will also be computed if nonzero */
	int do_gradient;

	/* periodic simulation box size */
	six_t box;

	/* stress tensor */
	mat_t stress;

	/* force and torque on fragments */
	six_t *grad;

	/* number of point charges or QM atoms */
	size_t n_ptc;

	/* coordinates of point charges or QM atoms */
	vec_t *ptc_xyz;

	/* point charges or QM atoms */
	double *ptc;

	/* gradient on point charges or QM atoms */
	vec_t *ptc_grad;

    /* total number of polarizable points */
	size_t n_polarizable_pts;

	/* number of core orbitals in ab initio subsystem */
	size_t n_ai_core;

	/* number of active orbitals in ab initio subsystem */
	size_t n_ai_act;

	/* number of virtual orbitals in ab initio subsystem */
	size_t n_ai_vir;

	/* ab initio orbital energies
	 * size [n_ai_occ + n_ai_vir] */
	double *ai_orbital_energies;

	/* ab initio dipole moment integrals on polarizable points
	 * size [3 * (n_ai_occ + n_ai_vir) ^ 2] */
	double *ai_dipole_integrals;

	/* EFP energy terms */
	struct efp_energy energy;

	/* skip-list of fragments - boolean array of nfrag^2 elements */
	char *skiplist;

	/* energies of pairwise interactions with the ligand */
	struct efp_energy *pair_energies;

	/* the number of symmetrically unique fragments. 0 for non-symmetric systems */
	size_t nsymm_frag;

	/* symmetry list - list of n_frag length specifying symmetrically identical fragments */
	size_t *symmlist;
};

#endif /* LIBEFP_PRIVATE_H */
