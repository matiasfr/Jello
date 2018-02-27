/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#include "jello.h"
#include "physics.h"

/* Computes acceleration to every control point of the jello cube,
   which is in state given by 'jello'.
   Returns result in array 'a'. */
void computeAcceleration(struct world * jello, struct point a[8][8][8])
{
  //Compute hooke's law //reminder F = ma ,  a = F/m
  //collision and non-collision k and damping constants
  double kE = jello->kElastic;
  double dE = jello->dElastic;
  double kC = jello->kCollision;
  double dC = jello->dCollision;

  point forces[8][8][8];

  point springForce = {}; //total force on starting point from given spring

//iterate through springs to calculate forces
  for(int i=0; i<springs.size(); i++) {
    //calculate hook force
    pSUM(springForce, calculateHooke(kE, springs[i].length, jello->p[springs[i].start[0]][springs[i].start[1]][springs[i].start[2]], jello->p[springs[i].end[0]][springs[i].end[1]][springs[i].end[2]]), springForce);

    //calculate damping force
    pSUM(springForce, calculateDamping(dE, jello->p[springs[i].start[0]][springs[i].start[1]][springs[i].start[2]], jello->p[springs[i].end[0]][springs[i].end[1]][springs[i].end[2]], jello->v[springs[i].start[0]][springs[i].start[1]][springs[i].start[2]], jello->v[springs[i].end[0]][springs[i].end[1]][springs[i].end[2]]), springForce);

    //add the accumulated forces to the point's force representation
    pSUM(forces[springs[i].start[0]][springs[i].start[1]][springs[i].start[2]], springForce, forces[springs[i].start[0]][springs[i].start[1]][springs[i].start[2]]);
  }


//use forces to calculate acceleration for each point
  for (int i=0; i<=7; i++) {
    for (int j=0; j<=7; j++) {
      for (int k=0; k<=7; k++) {
        pMULTIPLY(forces[springs[i].start[0]][springs[i].start[1]][springs[i].start[2]], (1.0 / jello->mass), a[i][j][k]);
      }
    }
  }


}

point calculateHooke(double k, double restLength, struct point a, struct point b) {
  point temp = {};
  temp.x = k * (-1.0) * (distanceHelper(a,b) - restLength) * ((a.x - b.x) / distanceHelper(a,b));
  temp.y = k * (-1.0) * (distanceHelper(a,b) - restLength) * ((a.y - b.y) / distanceHelper(a,b));
  temp.z = k * (-1.0) * (distanceHelper(a,b) - restLength) * ((a.z - b.z) / distanceHelper(a,b));
  return temp;
}

point calculateDamping(double k, struct point a, struct point b, struct point vA,struct point vB) {
  point temp = {};
  point velocityDiff = {};
  pDIFFERENCE(vA, vB, velocityDiff);
  double dampingCoefficient = ( velocityDiff.x * (a.x-b.x) + velocityDiff.y * (a.y-b.y) + velocityDiff.z * (a.z-b.z) ) / distanceHelper(a,b);

  temp.x = k * (-1.0) * dampingCoefficient * (a.x-b.x) / distanceHelper(a,b);
  temp.y = k * (-1.0) * dampingCoefficient * (a.y-b.y) / distanceHelper(a,b);
  temp.z = k * (-1.0) * dampingCoefficient * (a.z-b.z) / distanceHelper(a,b);

  return temp;
}


double distanceHelper(struct point a, struct point b) {
  return sqrt( pow((a.x-b.x), 2) + pow((a.y-b.y), 2) + pow((a.z-b.z), 2) );
}

void initializePoint(point &a) {
    a.x = 0;
    a.y = 0;
    a.z = 0;
}

/* performs one step of Euler Integration */
/* as a result, updates the jello structure */
void Euler(struct world * jello)
{
  int i,j,k;
  point a[8][8][8];

  computeAcceleration(jello, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
        jello->p[i][j][k].x += jello->dt * jello->v[i][j][k].x;
        jello->p[i][j][k].y += jello->dt * jello->v[i][j][k].y;
        jello->p[i][j][k].z += jello->dt * jello->v[i][j][k].z;
        jello->v[i][j][k].x += jello->dt * a[i][j][k].x;
        jello->v[i][j][k].y += jello->dt * a[i][j][k].y;
        jello->v[i][j][k].z += jello->dt * a[i][j][k].z;

      }
}

/* performs one step of RK4 Integration */
/* as a result, updates the jello structure */
void RK4(struct world * jello)
{
  point F1p[8][8][8], F1v[8][8][8],
        F2p[8][8][8], F2v[8][8][8],
        F3p[8][8][8], F3v[8][8][8],
        F4p[8][8][8], F4v[8][8][8];

  point a[8][8][8];


  struct world buffer;

  int i,j,k;

  buffer = *jello; // make a copy of jello

  computeAcceleration(jello, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         pMULTIPLY(jello->v[i][j][k],jello->dt,F1p[i][j][k]);
         pMULTIPLY(a[i][j][k],jello->dt,F1v[i][j][k]);
         pMULTIPLY(F1p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F1v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F2p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F2p[i][j][k]);
         // F2v = dt * a(buffer.p,buffer.v);
         pMULTIPLY(a[i][j][k],jello->dt,F2v[i][j][k]);
         pMULTIPLY(F2p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F2v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F3p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);
         pMULTIPLY(a[i][j][k],jello->dt,F3v[i][j][k]);
         pMULTIPLY(F3p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);


  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F4p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);
         pMULTIPLY(a[i][j][k],jello->dt,F4v[i][j][k]);

         pMULTIPLY(F2p[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3p[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1p[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4p[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->p[i][j][k],jello->p[i][j][k]);

         pMULTIPLY(F2v[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4v[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->v[i][j][k],jello->v[i][j][k]);
      }

  return;
}
