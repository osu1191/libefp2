/*-
 * Copyright (c) 2012-2015 Ilya Kaliman
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

#include "common.h"

#include <stream.h>

static void check_double(const struct cfg *cfg, const char *name)
{
	if (cfg_get_double(cfg, name) <= 0.0)
		error("%s must be greater than zero", name);
}

static void check_int(const struct cfg *cfg, const char *name)
{
	if (cfg_get_int(cfg, name) <= 0)
		error("%s must be greater than zero", name);
}

static void check_cfg(struct cfg *cfg)
{
	check_int(cfg, "max_steps");
	check_int(cfg, "print_step");
	check_double(cfg, "opt_tol");
	check_double(cfg, "num_step_dist");
	check_double(cfg, "num_step_angle");
	check_double(cfg, "time_step");
	check_double(cfg, "temperature");
	check_double(cfg, "pressure");
	check_double(cfg, "thermostat_tau");
	check_double(cfg, "barostat_tau");
}

static void parse_frag(struct stream *stream, enum efp_coord_type coord_type,
    struct frag *frag)
{
	memset(frag, 0, sizeof(struct frag));
	efp_stream_skip_space(stream);
	const char *ptr = efp_stream_get_ptr(stream);
	efp_stream_skip_nonspace(stream);
	size_t len = efp_stream_get_ptr(stream) - ptr;

	if (len == 0)
		error("fragment name is missing");
	frag->name = xmalloc(len + 1);

	for (size_t i = 0; i < len; i++)
		frag->name[i] = (char)tolower(*ptr++);
	frag->name[len] = '\0';

	efp_stream_next_line(stream);

	int n_rows = (int []) {
		[EFP_COORD_TYPE_XYZABC] = 1,
		[EFP_COORD_TYPE_POINTS] = 3,
		[EFP_COORD_TYPE_ROTMAT] = 4,
		[EFP_COORD_TYPE_ATOMS] = 3 }[coord_type];

	int n_cols = (int []) {
		[EFP_COORD_TYPE_XYZABC] = 6,
		[EFP_COORD_TYPE_POINTS] = 3,
		[EFP_COORD_TYPE_ROTMAT] = 3,
		[EFP_COORD_TYPE_ATOMS] = 4 }[coord_type];

	if (coord_type == EFP_COORD_TYPE_ATOMS) {
	    while (true) {
            struct efp_atom atom_i;
            memset(&atom_i, 0, sizeof(struct efp_atom));

            efp_stream_skip_space(stream);
            if (efp_strncasecmp(efp_stream_get_ptr(stream), "velocity", strlen("velocity")) == 0 ||
                efp_strncasecmp(efp_stream_get_ptr(stream), "constraint", strlen("constraint")) == 0 ||
                efp_strncasecmp(efp_stream_get_ptr(stream), "fragment", strlen("fragment")) == 0) {
                break;
            }

            const char *ptr2 = efp_stream_get_ptr(stream);
            efp_stream_skip_nonspace(stream);

            size_t len2 = efp_stream_get_ptr(stream) - ptr2;
            if (len2 == 0)
                error("incorrect fragment coordinates format: reading efp atom name");
            memcpy(atom_i.label, ptr2, len2);

            efp_stream_skip_space(stream);
            if (!efp_stream_parse_double(stream, &atom_i.x) || !efp_stream_parse_double(stream, &atom_i.y)
                || !efp_stream_parse_double(stream, &atom_i.z))
                error("incorrect fragment coordinates format: reading efp atom coordinates");	

            frag->n_atoms++;
            frag->atoms = xrealloc(frag->atoms, frag->n_atoms * sizeof(struct efp_atom));
            frag->atoms[frag->n_atoms - 1] = atom_i;
            
			frag->coord = xrealloc(frag->coord, frag->n_atoms * 3* sizeof(double));
            frag->coord[(frag->n_atoms-1)*3] = atom_i.x;
            frag->coord[(frag->n_atoms-1)*3+1] = atom_i.y;
            frag->coord[(frag->n_atoms-1)*3+2] = atom_i.z;

            //printf("n_atoms = %d",frag->n_atoms);

            efp_stream_next_line(stream);
            if (efp_stream_eol(stream))
                return;
        }
    }
	else {
		frag->coord = (double*)malloc(n_rows * n_cols * sizeof(double));
        for (int i = 0, idx = 0; i < n_rows; i++) {
            for (int j = 0; j < n_cols; j++, idx++) {
                if (!efp_stream_parse_double(stream, frag->coord + idx))
//		    printf("We are at line 140\n");
                    error("incorrect fragment coordinates format");
            }

            efp_stream_next_line(stream);
        }
    }

	efp_stream_skip_space(stream);
	if (efp_stream_eol(stream))
		return;

	if (efp_strncasecmp(efp_stream_get_ptr(stream), "velocity",
			strlen("velocity")) == 0) {
		efp_stream_next_line(stream);

		for (int i = 0; i < 6; i++)
			if (!efp_stream_parse_double(stream, frag->vel + i))
				error("incorrect fragment velocities format");

		efp_stream_next_line(stream);
	}

	efp_stream_skip_space(stream);
	if (efp_stream_eol(stream))
		return;

	if (efp_strncasecmp(efp_stream_get_ptr(stream), "constraint",
			strlen("constraint")) == 0) {
		efp_stream_next_line(stream);

		frag->constraint_enable = true;
		if (!efp_stream_parse_double(stream, &frag->constraint_k))
			error("incorrect fragment constraint format");
		if (!efp_stream_parse_double(stream, &frag->constraint_xyz.x))
			error("incorrect fragment constraint format");
		if (!efp_stream_parse_double(stream, &frag->constraint_xyz.y))
			error("incorrect fragment constraint format");
		if (!efp_stream_parse_double(stream, &frag->constraint_xyz.z))
			error("incorrect fragment constraint format");

		efp_stream_next_line(stream);
	}
}

static bool is_keyword(const char *str, const char *key)
{
	return efp_strncasecmp(str, key, strlen(key)) == 0 &&
	    isspace(str[strlen(key)]);
}

struct sys *parse_input(struct cfg *cfg, const char *path)
{
	struct stream *stream;
	struct sys *sys = xcalloc(1, sizeof(struct sys));

	if ((stream = efp_stream_open(path)) == NULL)
		error("unable to open input file");

	efp_stream_next_line(stream);

	while (!efp_stream_eof(stream)) {
		if (efp_stream_current_char(stream) == '#')
			goto next;

		efp_stream_skip_space(stream);

		if (efp_stream_eol(stream))
			goto next;

		if (is_keyword(efp_stream_get_ptr(stream), "fragment")) {
			struct frag frag;
			enum efp_coord_type coord_type;

			efp_stream_advance(stream, strlen("fragment"));
			coord_type = cfg_get_enum(cfg, "coord");
			parse_frag(stream, coord_type, &frag);

			sys->n_frags++;
			sys->frags = xrealloc(sys->frags,
			    sys->n_frags * sizeof(struct frag));
			sys->frags[sys->n_frags - 1] = frag;
			continue;
		}
		else if (is_keyword(efp_stream_get_ptr(stream), "charge")) {
			struct charge charge;

			memset(&charge, 0, sizeof(charge));
			efp_stream_advance(stream, strlen("charge"));

			if (!efp_stream_parse_double(stream, &charge.q) ||
			    !efp_stream_parse_double(stream, &charge.pos.x) ||
			    !efp_stream_parse_double(stream, &charge.pos.y) ||
			    !efp_stream_parse_double(stream, &charge.pos.z))
				error("unable to parse a point charge record");

			sys->n_charges++;
			sys->charges = xrealloc(sys->charges,
			    sys->n_charges * sizeof(struct charge));
			sys->charges[sys->n_charges - 1] = charge;
		}
		else {
			if (!cfg_parse_line(cfg, efp_stream_get_ptr(stream)))
				error("%s", cfg_get_last_error(cfg));
		}
next:
		efp_stream_next_line(stream);
	}

	if (sys->n_frags == 0)
		error("at least one fragment must be specified");

	if (sys->n_charges > 0 && cfg_get_enum(cfg, "run_type") == RUN_TYPE_MD)
		error("point charges are not supported in molecular dynamics");

	if (cfg_get_enum(cfg, "ensemble") == ENSEMBLE_TYPE_NPT)
		cfg_set_bool(cfg, "enable_pbc", true);

	if (cfg_get_bool(cfg, "enable_pbc"))
		cfg_set_bool(cfg, "enable_cutoff", true);

	// turn off pairwise calculations if ligand is not set
	if (cfg_get_int(cfg, "ligand") == -100 && cfg_get_bool(cfg, "enable_pairwise") == true)
	    error("Specify ligand for pairwise calculations");

	check_cfg(cfg);
	efp_stream_close(stream);

	return sys;
}
