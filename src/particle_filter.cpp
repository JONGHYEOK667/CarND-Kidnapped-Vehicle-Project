/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;
using std::uniform_real_distribution;
using std::uniform_int_distribution;


void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 200;  // TODO: Set the number of particles
   
  double gps_x = x;
  double gps_y = y;
  double gps_theta = theta;

  double std_x = std[0];
  double std_y = std[1];
  double std_theta = std[2];


  std::cout << "--Particle Filter Init info--" << "\n";
  std::cout << "num_particles : " << num_particles << "\n";
  std::cout << "gps_x : " << gps_x << "\n";
  std::cout << "gps_y : " << gps_y << "\n";
  std::cout << "gps_theta : " << gps_theta << "\n";
  std::cout << "std_x : " << std_x << "\n";
  std::cout << "std_y : " << std_y << "\n";
  std::cout << "std_theta : " << std_theta << "\n";

  
 
  normal_distribution<double> dist_x(gps_x, std_x);
  normal_distribution<double> dist_y(gps_y, std_y);
  normal_distribution<double> dist_theta(gps_theta, std_theta);

  std::default_random_engine gen;

  for (int i = 0; i < num_particles; ++i) {
      Particle particle;

      particle.id = i;
      particle.x = dist_x(gen);
      particle.y = dist_y(gen);
      particle.theta = dist_theta(gen);
      particle.weight = 1.0;


      
      if (i < 2) {
          std::cout << "--first 2 particle Init info--" << "\n";
          std::cout << "particle ID : " << particle.id <<"\n";
          std::cout << "particle x : " << particle.x << "\n";
          std::cout << "particle y : " << particle.y << "\n";
          std::cout << "particle theta : " << particle.theta << "\n";
      }

      particles.push_back(particle);

  }

  is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */

  double std_x = std_pos[0];
  double std_y = std_pos[1];
  double std_theta = std_pos[2];

  normal_distribution<double> dist_x(0, std_x);
  normal_distribution<double> dist_y(0, std_y);
  normal_distribution<double> dist_theta(0, std_theta);

  std::default_random_engine gen;

  /*
  std::cout << "--Velocity and Yawrate for Predict--" << "\n";
  std::cout << "Velocity : " << velocity << "\n";
  std::cout << "yaw_rate : " << yaw_rate << "\n";
  */



  for (int i = 0; i < num_particles; ++i) {


      if (fabs(yaw_rate) < 0.00001) {

          particles[i].x += velocity * delta_t * cos(particles[i].theta);;
          particles[i].y += velocity * delta_t * cos(particles[i].theta);;
      
      }
      else {

          particles[i].x = particles[i].x + velocity / yaw_rate * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta));
          particles[i].y = particles[i].y + velocity / yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
          particles[i].theta += yaw_rate * delta_t;

      }
      



      /*
      std::cout << "--first 2 particle Predict info--" << "\n";
      if (i < 2) {
          std::cout << "particle ID : " << particles[i].id << "\n";
          std::cout << "particle x : " << particles[i].x << "\n";
          std::cout << "particle y : " << particles[i].y << "\n";
          std::cout << "particle theta : " << particles[i].theta << "\n";
      }
      */

      particles[i].x += dist_x(gen);
      particles[i].y += dist_y(gen);
      particles[i].theta += dist_theta(gen);

      

      /*
      std::cout << "--first 2 particle Predict Add noise info--" << "\n";
      if (i < 2) {
          std::cout << "particle ID : " << particles[i].id << "\n";
          std::cout << "particle x : " << particles[i].x << "\n";
          std::cout << "particle y : " << particles[i].y << "\n";
          std::cout << "particle theta : " << particles[i].theta << "\n";
      }
      */


  }

  


}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */

  int nPredicted = predicted.size();
  int nObservations = observations.size();

  for (int i = 0; i < nObservations; ++i) {

      double minRange = std::numeric_limits<double>::max();
      int Lm_id = -1;

      for (int j = 0; j < nPredicted; ++j) {
      
          double dx = observations[i].x - predicted[j].x;
          double dy = observations[i].y - predicted[j].y;

          double dr = dx * dx + dy * dy;

          if (dr < minRange) {

              minRange = dr;
              Lm_id = predicted[j].id;

          }      
      }

      observations[i].id = Lm_id; 
  
  }
  

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

  double std_x = std_landmark[0];
  double std_y = std_landmark[1];

  for (int i = 0; i < num_particles; ++i) {
      
      double x = particles[i].x;
      double y = particles[i].y;
      double theta = particles[i].theta;

      /*
      std::cout << "i : " << i << std::endl;
      std::cout << "x : " << x << std::endl;
      std::cout << "y : " << y << std::endl;
      std::cout << "theta : " << theta << std::endl;
      */


      vector<LandmarkObs> inRangeLandmarks;
      for (int j = 0; j < map_landmarks.landmark_list.size(); ++j) {

          double x_Lm = map_landmarks.landmark_list[j].x_f;
          double y_Lm = map_landmarks.landmark_list[j].y_f;
          double id_Lm = map_landmarks.landmark_list[j].id_i;

          double delta_x = x - x_Lm;
          double delta_y = y - y_Lm;

          if ((sensor_range * sensor_range) >= (delta_x * delta_x) + (delta_y * delta_y)) {
              /*
              std::cout << "x_Lm : " << x_Lm << std::endl;
              std::cout << "y_Lm : " << y_Lm << std::endl;
              std::cout << "id_Lm : " << id_Lm << std::endl;
              */
              inRangeLandmarks.push_back(LandmarkObs{ id_Lm, x_Lm, y_Lm });

          }
      }
      vector<LandmarkObs> observedLandmarks;
      for (int j = 0; j < observations.size(); ++j) {
          
          double x_Ob = observations[j].x;
          double y_Ob = observations[j].y;
          int id_Ob = observations[j].id;

          double x_Ob_Map = cos(theta) * x_Ob - sin(theta) * y_Ob + x;
          double y_Ob_Map = sin(theta) * x_Ob + cos(theta) * y_Ob + y;
          
          /*
          std::cout << "x_Ob : " << x_Ob << std::endl;
          std::cout << "y_Ob : " << y_Ob << std::endl;
          std::cout << "id_Ob : " << id_Ob << std::endl;
          */
          observedLandmarks.push_back(LandmarkObs{ id_Ob, x_Ob_Map, y_Ob_Map });

      }

      dataAssociation(inRangeLandmarks, observedLandmarks);

      particles[i].weight = 1.0;
      for (int j = 0; j < observedLandmarks.size(); j++) {

          double x_Ob_Map = observedLandmarks[j].x;
          double y_Ob_Map = observedLandmarks[j].y;
          int id_Ob_Map = observedLandmarks[j].id;


          double x_Lm, y_Lm;
          int id_Lm;
          int k = 0;
          int nLandmark = inRangeLandmarks.size();
          bool isFound = false;

          while (!isFound && k < nLandmark) {

              if (inRangeLandmarks[k].id == id_Ob_Map) {

                  isFound = true;
                  x_Lm = inRangeLandmarks[k].x;
                  y_Lm = inRangeLandmarks[k].y;
                  id_Lm = inRangeLandmarks[k].id;
                               
              }

              k += 1;

          }
          /*
          std::cout << "j : " << j << std::endl;
          std::cout << "x_Ob_Map : " << x_Ob_Map << std::endl;
          std::cout << "y_Ob_Map : " << y_Ob_Map << std::endl;
          std::cout << "id_Ob_Map : " << id_Ob_Map << std::endl;
          */
          /*
          std::cout << "x_Lm : " << x_Lm << std::endl;
          std::cout << "y_Lm : " << y_Lm << std::endl;
          std::cout << "id_Lm : " << id_Lm << std::endl;
          */
          double dx = x_Ob_Map - x_Lm;
          double dy = y_Ob_Map - y_Lm;
          double weight_temp;
          weight_temp = (1 / (2 * M_PI * std_x * std_y)) * exp(-(dx * dx / (2 * std_x * std_x) + (dy * dy / (2 * std_y * std_y))));
          particles[i].weight *= weight_temp;

          
      }
  }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

  vector<double> weights;
  double maxWeight = std::numeric_limits<double>::min();
  std::default_random_engine gen;

  for (int i = 0; i < num_particles; i++) {

      weights.push_back(particles[i].weight);

      if (particles[i].weight > maxWeight) {
          maxWeight = particles[i].weight;
      }

  }

  uniform_real_distribution<double> distDouble(0.0, maxWeight);
  uniform_int_distribution<int> distInt(0, num_particles - 1);
  int index = distInt(gen);
  double beta = 0.0;

  vector<Particle> resampledParticles;
  for (int i = 0; i < num_particles; i++) {
      beta += distDouble(gen) * 2.0;
      while (beta > weights[index]) {
          beta -= weights[index];
          index = (index + 1) % num_particles;
      }
      resampledParticles.push_back(particles[index]);
  }

  particles = resampledParticles;

}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}