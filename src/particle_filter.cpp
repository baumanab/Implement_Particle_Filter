/*
 * particle_filter.cpp
 *
 *  Template Code Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include <limits> // support min/max searching
#include <cmath> // supports math....obviously

#include "particle_filter.h"

using namespace std;

// initialize random number generator to be used globally
default_random_engine gen;


void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	// x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    // NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    // NOTE: This ties directly to this quiz imp part filter concepts 4, 5  	


    num_particles = 100;
    
	// generate normal distribution with mean 0 std_dev from indexed std vector
	// std vector is entered in main.cpp as std_pos
    normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);
	
	
    // cycle each particle, initialize by sampling, add noise
	for (int i = 0; i < num_particles; i++) {
		Particle particle;
        
		// set attributes and sample from distribution
		particle.id = i;
		particle.x = dist_x(gen);
		particle.y = dist_y(gen);
		particle.theta = dist_theta(gen);
		particle.weight = 1.0;
        
		// push_back adds element to end of vector and expands length by 1
		// https://stackoverflow.com/questions/13324431/c-vectors-insert-push-back-difference
		particles.push_back(particle);
		weights.push_back(particle.weight);
	}
	// switch bool
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	// http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	// http://www.cplusplus.com/reference/random/default_random_engine/	

    // cycle each particle, calculate new positons given +/- yaw, add noise
	for (int i = 0; i < num_particles; i++) {

		double new_x, new_y, new_theta;
		
		if((yaw_rate < 0.0001) & yaw_rate > -0.0001) { // 0 +/- 0.001, so 0 within error
		    new_theta = particles[i].theta;
			new_x = particles[i].x + (velocity * delta_t * cos(new_theta));
			new_y = particles[i].y + (velocity * delta_t * sin(new_theta));
		
		}
		else {	
			new_theta = particles[i].theta + yaw_rate * delta_t;
			new_x = particles[i].x + ((velocity/yaw_rate) * (sin(new_theta) - sin(particles[i].theta)));
			new_y = particles[i].y + ((velocity/yaw_rate) * (cos(particles[i].theta) - cos(new_theta)));

		}

		normal_distribution<double> dist_x(new_x, std_pos[0]);
		normal_distribution<double> dist_y(new_y, std_pos[1]);
		normal_distribution<double> dist_theta(new_theta, std_pos[2]);

		// sample from distribution to add noise
		particles[i].x = dist_x(gen);
		particles[i].y = dist_y(gen);
		particles[i].theta = dist_theta(gen);

	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
    // observed measurement to this particular landmark.
    // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
    // implement this method and use it as a helper during the updateWeights phase.
    // NOTE: we have predicted measurement vector and an observations vector
    // iterate over observations calculating distance between predicted and observed
    // take min of either by track and update or creating diff vector and extrating min value
    // NOTE:  This calls observations by reference (>&) so this function modifies the original vector

	for (int i = 0; i < observations.size(); i++) {
        double min_dist = numeric_limits<double>::infinity(); // largest possible value
		for (int j = 0; j < predicted.size(); j++) {
			double current_dist_ = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
			if (current_dist_ < min_dist) {
				min_dist = current_dist_;
				observations[i].id = predicted[j].id;
			}
		}
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],std::vector<LandmarkObs>& observations, Map& map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
    // more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
    // NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
    // according to the MAP'S coordinate system. You will need to transform between the two systems.
    // Keep in mind that this transformation requires both rotation AND translation (but no scaling).
    // The following is a good resource for the theory:
    // https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
    // and the following is a good resource for the actual equation to implement (look at equation
    // 3.33 http://planning.cs.uiuc.edu/node99.html
    // NOTE: We have an landmark observations vector in vehicle coords, map_landmarks,
    // sensor_range, and stdev for the landmark.
    // make conversion, make association using the func above --> mv-gauss
    // mulitply result to get final weight
    // NOTE:  What do I need: obervations[x, y], particle[x, y, theta]
	
	// extract map landmark (x,y) and observations (x,y) and stdev
    double sig_x = std_landmark[0];
    double sig_y = std_landmark[1];
	double gauss_norm = 2 * M_PI * sig_x * sig_y;    
	double total_weight = 0.0;
	
    for (int i = 0; i < num_particles; i++) {

        double particle_x = particles[i].x;
        double particle_y = particles[i].y;
        double particle_theta = particles[i].theta;
		
		// initialize a vector as container for transformed coordinate pairs
		vector<LandmarkObs> obs_map_coords;
		
		// homogenous transform to map coordinates
		for (int j = 0; j < observations.size(); j++) {
			LandmarkObs obs_map;
			double obs_x = observations[j].x;
			double obs_y = observations[j].y;
			// using anti-clockwise rotation matrix
			obs_map.x = particle_x + (obs_x * cos(particle_theta) - obs_y * sin(particle_theta));
			obs_map.y = particle_y + (obs_x * sin(particle_theta) + obs_y * cos(particle_theta));
			obs_map_coords.push_back(obs_map);
		}
		
		// detectable landmarks == landmarks in range of sensor_range
        // container for predicted landmarks
		vector<LandmarkObs> detectable_landmarks;
		
		// search for landmarks within sensor range
        for (int k = 0; k < map_landmarks.landmark_list.size(); k++) {
			double dist_ = dist(particle_x, particle_y, map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f);
			
			if (dist_ < sensor_range) {
				LandmarkObs detected_landmarks;
				detected_landmarks.id = map_landmarks.landmark_list[k].id_i;
				detected_landmarks.x = map_landmarks.landmark_list[k].x_f;
				detected_landmarks.y = map_landmarks.landmark_list[k].y_f;

				detectable_landmarks.push_back(detected_landmarks);
			}
		}

        // associate observation in map coords with closest predicted landmark
        // NOTE: In this function the observations vector is called by reference
        // so this function modifies the original vector, for use in search

        // associate map_landmarks id and observations in map coordinates id to be used as foreign key

		dataAssociation(detectable_landmarks, obs_map_coords);
		
        double updated_weight = 1.0;
        double mu_x, mu_y; // mu is actual landmark location (ground truth)
		particles[i].weight = 1.0;
		
		// extract map landmark (x,y) and observations (x,y) and stdev
		for (int j = 0; j < obs_map_coords.size(); j++) {

			for (int k = 0; k < detectable_landmarks.size(); k++) {

				// find the actual landmark that corresponds to the observation map coords
				if (obs_map_coords[j].id == detectable_landmarks[k].id) {
					mu_x = detectable_landmarks[k].x;
					mu_y = detectable_landmarks[k].y;
					break;
				}
			}
			
			double x_obs = obs_map_coords[j].x;
            double y_obs = obs_map_coords[j].y;
			double x_diff = x_obs - mu_x;
            double y_diff = y_obs - mu_y;
			double x_delta = pow((x_diff)/(2 * sig_x), 2.0);
			double y_delta = pow((y_diff)/(2 * sig_y), 2.0);

			double gauss_prob = exp(-(x_delta + y_delta)) * (1/gauss_norm);
			updated_weight *= gauss_prob;			

		}

        particles[i].weight = updated_weight;
		weights[i] = updated_weight;
		total_weight += updated_weight;
    }

	// normalize weight vector
	if (total_weight != 0) { //let's not divide by 0 :)
		for (int i = 0; i < weights.size(); i++) {
			particles[i].weight /= total_weight;
			weights[i] /= total_weight;
		}
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	// http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    // this solution is from walk-through.  Could also implement a sampling wheel
	
	discrete_distribution<int> distrib(weights.begin(), weights.end());
	vector<Particle> resampled_particles;
	
	for (int i = 0; i < num_particles; i++) {
		resampled_particles.push_back(particles[distrib(gen)]);
	}
	particles = resampled_particles;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
