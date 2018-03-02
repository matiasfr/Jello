/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#ifndef _PHYSICS_H_
#define _PHYSICS_H_

void computeAcceleration(struct world * jello, struct point a[8][8][8]);
void initializePoint(point &a);

// perform one step of Euler and Runge-Kutta-4th-order integrators
// updates the jello structure accordingly
void Euler(struct world * jello);
void RK4(struct world * jello);
double distanceHelper(struct point a, struct point b);
point calculateHooke(double k, double restLength, struct point a, struct point b);
point calculateDamping(double k, struct point a, struct point b, struct point vA,struct point vB);
void normalize(struct point & diff);
#endif
