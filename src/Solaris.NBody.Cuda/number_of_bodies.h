#pragma once

#include "interaction_bound.h"

class number_of_bodies {
public:
	number_of_bodies(int star, int giant_planet, int rocky_planet, int proto_planet, int super_planetesimal, int planetesimal, int test_particle);

	//! Calculates the number of bodies with mass, i.e. sum of the number of stars, giant planets, 
	/*  rocky planets, protoplanets, super-planetesimals and planetesimals.
	*/
	int		n_massive();
	//! Calculates the number of bodies which are self-interacting, i.e. sum of the number of stars, giant planets, 
	/*  rocky planets and protoplanets.
	*/
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
