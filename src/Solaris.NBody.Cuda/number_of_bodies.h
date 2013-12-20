#pragma once

#include "interaction_bound.h"

class number_of_bodies {
public:
	number_of_bodies(int star, int giant_planet, int rocky_planet, int proto_planet, int super_planetesimal, int planetesimal, int test_particle);

	int		n_massive();
	int		n_self_interacting();

	interaction_bound get_self_interacting();
	interaction_bound get_nonself_interacting();
	interaction_bound get_non_interacting();
	interaction_bound get_bodies_gasdrag();

	int		star;
	int		giant_planet;
	int		rocky_planet;
	int		proto_planet;
	int		super_planetesimal;
	int		planetesimal;
	int		test_particle;
	int		total;
	int		total_rounded_up;

private:
	int2_t		sink;
	int2_t		source;
};
