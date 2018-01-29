/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
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
#include <vector>
#include <limits>

#include "particle_filter.h"

using namespace std;
static default_random_engine gen;
void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 101;

	normal_distribution<double> dist_xi(0,std[0]);
	normal_distribution<double> dist_yi(0,std[1]);
	normal_distribution<double> dist_thetai(0,std[2]);

	for(int i=0;i<num_particles;i++){
		Particle p;
		p.id = i;
		p.x = x + dist_xi(gen);
		p.y = y + dist_yi(gen);
		p.theta = theta + dist_thetai(gen);	
		p.weight = 1.0;
		particles.push_back(p);
	}
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	normal_distribution<double> dist_x(0,std_pos[0]);
	normal_distribution<double> dist_y(0,std_pos[1]);
	normal_distribution<double> dist_theta(0,std_pos[2]);

	for(int i=0;i<num_particles;i++){
		if(fabs(yaw_rate)<0.00001){
			particles[i].x = particles[i].x + velocity * delta_t * cos(particles[i].theta) + dist_x(gen);
			particles[i].y = particles[i].y + velocity * delta_t * sin(particles[i].theta) + dist_y(gen);
			particles[i].theta = particles[i].theta + dist_theta(gen);		
		}
		else{
			particles[i].x = particles[i].x + (velocity/yaw_rate)*(sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta)) + dist_x(gen);
			particles[i].y = particles[i].y + (velocity/yaw_rate)*(-cos(particles[i].theta + yaw_rate*delta_t) + cos(particles[i].theta)) + dist_y(gen);
			particles[i].theta = particles[i].theta + yaw_rate*delta_t + dist_theta(gen);	
		}	
	}
	
}

//void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
void ParticleFilter::dataAssociation(std::vector<LandmarkObs> map_landmarks_in_range, std::vector<LandmarkObs>& observations_in_map_coordinates) {
	// TODO: Find the predicted measurement that is closest to each observed measurement (map) and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	std::vector<LandmarkObs> result;

	for (int i=0;i<observations_in_map_coordinates.size();i++){
		double dist1 = std::numeric_limits<double>::max();
		int map_index = -1;
		double x_up,y_up,minimum;

		for (int j = 0;j<map_landmarks_in_range.size();j++){
			minimum = dist(map_landmarks_in_range[j].x, map_landmarks_in_range[j].y, observations_in_map_coordinates[i].x, observations_in_map_coordinates[i].y);
			if(minimum < dist1){
				dist1 = minimum;
				map_index = map_landmarks_in_range[j].id;
				//x_up = map_landmarks_in_range[j].x;
				//y_up = map_landmarks_in_range[j].y;
			}
		}
		//LandmarkObs l = {map_index,x_up,y_up};
		//result.push_back(l);
		observations_in_map_coordinates[i].id = map_index;
	}
	//observations_in_map_coordinates = result;
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	for (int i=0;i<num_particles;i++){
		double px = particles[i].x;
		double py = particles[i].y;
		double ptheta = particles[i].theta;
		
		//Finding the observations in range of sensor
		std::vector<LandmarkObs> map_landmarks_in_range;
		int l_id;
		double lx,ly;
		for(int j=0;j<map_landmarks.landmark_list.size();j++){
		l_id = map_landmarks.landmark_list[j].id_i;
		lx = map_landmarks.landmark_list[j].x_f;
		ly = map_landmarks.landmark_list[j].y_f;
			if(dist(px,py,lx,ly) <= sensor_range){
				LandmarkObs l = {l_id,lx,ly};
				map_landmarks_in_range.push_back(l);			
			}
		}
		//Transform observations to map coordinates
		std::vector<LandmarkObs> observations_in_map_coordinates;
		double xm,ym,ox,oy;
		int o_id;
		for (int j=0;j<observations.size();j++){
			o_id = observations[j].id;
			ox = observations[j].x;
			oy = observations[j].y;
			xm = ox*cos(ptheta) - oy*sin(ptheta) + px;
			ym = ox*sin(ptheta) + oy*cos(ptheta) + py;	
			LandmarkObs l = {o_id,xm,ym};
			observations_in_map_coordinates.push_back(l);
		}
		//Data Association
		dataAssociation(map_landmarks_in_range, observations_in_map_coordinates);	
		
		particles[i].weight = 1.0;
		double mvgpdf,mu_x,mu_y,o_x,o_y,nearest_landmark_index;
		for (int j=0;j<observations_in_map_coordinates.size();j++){
			o_x = observations_in_map_coordinates[j].x - px;
			o_y = observations_in_map_coordinates[j].y - py;
			nearest_landmark_index = observations_in_map_coordinates[j].id;
			
			for (int k = 0; k < map_landmarks_in_range.size(); k++){
        			if (map_landmarks_in_range[k].id == nearest_landmark_index){
          				mu_x = map_landmarks_in_range[k].x - px;
          				mu_y = map_landmarks_in_range[k].y - py;
        			}
			}
			//Multi - variate Gaussian Probability Density Function
			mvgpdf = (1/(2*M_PI*std_landmark[0]*std_landmark[1]))*exp(-(pow(((o_x-mu_x)/2*pow(std_landmark[0],2)),2)+pow(((o_y-mu_y)/2*pow(std_landmark[1],2)),2)));
			particles[i].weight = particles[i].weight*mvgpdf;
		}
	}
	double total_weight = 0;
	for (int i=0;i<num_particles;i++){
		total_weight = total_weight + particles[i].weight;
	}
	for (int i=0;i<num_particles;i++){
		particles[i].weight = particles[i].weight/total_weight;
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	double maximum=0.0;
	vector<Particle> particles_new;
	for (int i=0;i<num_particles;i++){
		if(particles[i].weight>maximum){
			maximum = particles[i].weight;		
		}
	} 
	srand((unsigned)time(0)); 
	//int index = int(((float)rand()/(float)RAND_MAX)*num_particles);
	uniform_int_distribution<int> index_calc(0,num_particles-1);
	int index = index_calc(gen);
	
	uniform_real_distribution<double> ramdom_beta(0.0, 1.0);
	double beta = 0.0;
	for (int i=0;i<num_particles;i++){
		beta = beta + ramdom_beta(gen)*2.0*maximum;
		while(particles[index].weight<beta){
			beta = beta - particles[index].weight;
			index = (index + 1)%num_particles;
		}
		particles_new.push_back(particles[index]);	
	}
	particles = particles_new; 
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
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
