/*
 * particle_filter.cpp
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

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	// x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	
	if(is_initialized) {
		return;
	} 
	num_particles = 50;

	default_random_engine rand_gen; 

	//create the normal distribution of the particles.
	normal_distribution<double> dist_x(x,std[0]);
	normal_distribution<double> dist_y(y,std[1]);
	normal_distribution<double> dist_t(theta,std[2]);

	// create particles and set the weights to 1
	// http://www.cplusplus.com/reference/random/default_random_engine/
	
	
	for(int i = 0 ; i < num_particles;i++){
		Particle p; //= new Particle();
		p.id = i;
		p.x = dist_x(rand_gen);
		p.y = dist_y(rand_gen);
		p.theta = dist_t(rand_gen);

		p.weight = 1.0;
		particles.push_back(p);
	}


	is_initialized = true;

	//cout << "\t init Done!" << endl;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	// http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	// http://www.cplusplus.com/reference/random/default_random_engine/
	
	//create the normal distribution for the sensor.
	//cout << "@ prediction starts!" << endl;

	normal_distribution<double> dist_x(0,std_pos[0]);
	normal_distribution<double> dist_y(0,std_pos[1]);
	normal_distribution<double> dist_t(0,std_pos[2]);

	default_random_engine rand_gen; 
	for(int i = 0 ; i < num_particles; i++){
		if( fabs(yaw_rate) < EPS){
			particles[i].x += velocity * delta_t * cos(particles[i].theta);
			particles[i].y += velocity * delta_t * sin(particles[i].theta);

		}
		else {

			particles[i].x += velocity/yaw_rate * (sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
			particles[i].y += velocity / yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
			particles[i].theta += yaw_rate * delta_t;
		}

		// add noise 
		particles[i].x += dist_x(rand_gen);
		particles[i].y += dist_y(rand_gen);
		particles[i].theta += dist_t(rand_gen);	
	}
	//cout << "\tprediction ends!" << endl;
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	//cout << "@ dataAssociation starts!" << endl;
	
	for(unsigned int i = 0 ; i < observations.size();i++){
		int id = -1;
		double min_dist = numeric_limits<double>::max();
		for(unsigned int j = 0 ; j < predicted.size();j++){
			double dist_cal = dist(observations[i].x,observations[i].y,predicted[j].x,predicted[j].y); //from helper functions
			if(dist_cal<min_dist){
				min_dist = dist_cal;
				id = predicted[j].id;
			}
		}
		observations[i].id = id;
	}
	//cout << "\tdataAssociation ends!" << endl;
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

	//cout << "@ updateWeights starts!" << endl;
	for(int i = 0 ; i < num_particles; i++){

		vector<LandmarkObs> near_landmarks ; 
		//update for each landmark
		for(unsigned int j = 0 ; j <  map_landmarks.landmark_list.size();j++){
			float landmark_x = map_landmarks.landmark_list[j].x_f;
			float landmark_y = map_landmarks.landmark_list[j].y_f;
			int landmark_id = map_landmarks.landmark_list[j].id_i;
			
			double dist_ = dist(particles[i].x,particles[i].y,
								landmark_x, landmark_y);

			if(dist_ <= sensor_range){
				near_landmarks.push_back(LandmarkObs{landmark_id,landmark_x,landmark_y}); //LandmarkObs defined in helper.h
			}	
		}

		//now transform the cordinates
		vector<LandmarkObs> translated_observations;
	    for(unsigned int j = 0; j < observations.size(); j++) {
	      double t_x = cos(particles[i].theta)*observations[j].x - sin(particles[i].theta)*observations[j].y + particles[i].x;
	      double t_y = sin(particles[i].theta)*observations[j].x + cos(particles[i].theta)*observations[j].y + particles[i].y;
	      translated_observations.push_back(LandmarkObs{observations[j].id, t_x, t_y});
	    }

	    //perform data association
	    dataAssociation(near_landmarks, translated_observations);

	    //reassign weight
	    particles[i].weight = 1.0;


	    // calculate weight
	    for(unsigned int j = 0 ; j < translated_observations.size(); j++){
	    	double o_x, o_y, pr_x, pr_y;
			o_x = translated_observations[j].x;
			o_y = translated_observations[j].y;
			int asso_prediction_id = translated_observations[j].id;

			//x,y coordinates of the prediction associated with the current observation
			for(unsigned int k = 0; k < near_landmarks.size(); k++) {
				if(near_landmarks[k].id == asso_prediction_id) {
					pr_x = near_landmarks[k].x;
					pr_y = near_landmarks[k].y;
				}
			}

			//Weight for this observation with multivariate Gaussian
			double s_x = std_landmark[0];
			double s_y = std_landmark[1];
			double obs_w = ( 1/(2*M_PI*s_x*s_y)) * exp( -( pow(pr_x-o_x,2)/(2*pow(s_x, 2)) + (pow(pr_y-o_y,2)/(2*pow(s_y, 2))) ) );

			//Product of this obersvation weight with total observations weight
			particles[i].weight *= obs_w;
	    }

	}
	//cout << "\tupdateWeights ends!" << endl;
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	//cout << "@ resample starts!" << endl;
	vector<double> weights; 
	double max_weight = numeric_limits<double>::min();
	for(int i = 0; i < num_particles; i++) {
		weights.push_back(particles[i].weight);
		if(particles[i].weight > max_weight) {
			max_weight = particles[i].weight;
		}
	}

	default_random_engine rand_gen; 
	uniform_real_distribution<double> distDouble(0.0, max_weight);
	uniform_int_distribution<int> distInt(0, num_particles - 1);


	int index = distInt(rand_gen);
	double beta = 0.0;
	vector<Particle> resampled_particles;
	for(int i = 0; i < num_particles; i++) {
		beta += distDouble(rand_gen) * 2.0;
		while(beta > weights[index]) {
			beta -= weights[index];
			index = (index + 1) % num_particles;
		}
		resampled_particles.push_back(particles[index]);
	}

	particles = resampled_particles;

	//cout << "\tresample ends!" << endl;

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
