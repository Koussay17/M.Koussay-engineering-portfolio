/** Title: Taylor-Culick retraction with relaxation...
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
# Updated: Oct 7 2023
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#define FILTERED
#include "two-phase.h"

#define ELASTIC 1

#if ELASTIC
// here, \sigma_p at the interface is ZERO
// #include "log-conform-PurelyElastic_v0.h"
// here, \sigma_p at the interface is NON-ZERO -> it comes from the stress balance from the momentum equation done at the interface
#include "log-conform-PurelyElastic_v1.h"
#else
// here, \sigma_p at the interface is ZERO
// #include "log-conform-ViscoElastic_v0.h"
// here, \sigma_p at the interface is NON-ZERO -> it comes from the stress balance from the momentum equation done at the interface
#include "log-conform-ViscoElastic_v1.h"
#endif

#include "navier-stokes/conserving.h"
#include "tension.h"

#define tsnap (0.25)
#define hole0 (1e0)
#define h0 (1e0)

// error tolerances
#define fErr (1e-3)
#define VelErr (1e-6)
#define KErr (1e-6)
#define AErr (1e-4)

double De, Ohf, Ohs, Ec, RhoR, tmax, Ldomain;
int MAXlevel;
scalar Gpd[], lambdav[];

int main(int argc, char const *argv[]) {
  dtmax = 1e-5;
  Ohs = 1e-5;
  RhoR = 1e-3;

  Ohf = 0.05;
  
  Ec = atof(argv[1]);
  // De = atof(argv[2]);
  MAXlevel = atoi(argv[3]);

  tmax = 25.0;
  Ldomain = 100.0;

  L0=Ldomain;
  X0=0.0; Y0=0.;
  init_grid (1 << (10));

  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);

  rho1 = 1.000; mu1 = Ohf;
  rho2 = RhoR; mu2 = Ohs;

  // polymers
  Gp = Gpd;
  // lambda = lambdav;

  f.sigma = 1.0;

  fprintf(ferr, "Capillary Scaling: Level %d tmax %g. Ohf %3.2e, Ec %3.2e, De Infty\n", MAXlevel, tmax, Ohf, Ec);
  run();
}

event properties (i++) {
  foreach () {
    Gpd[] = clamp(f[], 0., 1.)*Ec;
    // lambdav[] = clamp(f[], 0., 1.)*De;
  }
}

event init(t = 0){
  if(!restore (file = "dump")){
    /**
    We can now initialize the volume fractions in the domain. */
    refine(x<(h0/2.0)+0.025 && y < hole0+(h0/2.0)+0.025 && level<MAXlevel);
    fraction(f, y < hole0+(h0/2.0) ? sq(h0/2.0)-(sq(x)+sq(y-(h0/2.0)-hole0)) : (h0/2.0)-x);
  }
  // mask (x > 10.0 ? right : none);
}

scalar KAPPA[];
event adapt(i++) {
  curvature(f, KAPPA);
  // boundary((scalar *){KAPPA});
  adapt_wavelet ((scalar *){f, u.x, u.y, KAPPA, conform_p.x.x, conform_p.x.y, conform_p.y.y, conform_qq},
    (double[]){fErr, VelErr, VelErr, KErr, AErr, AErr, AErr, AErr},
    MAXlevel);
  unrefine(x>150.0);
}

// Outputs
event writingFiles (t = 0; t += tsnap; t <= tmax + tsnap) {
  dump (file = "dump");
  char nameOut[80];
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump (file = nameOut);
}

event logWriting (i++) {
  
  double ke = 0.;
  foreach (reduction(+:ke)){
    ke += sq(Delta)*(sq(u.y[])+sq(u.x[]))*f[];
  }

  static FILE * fp;
  if (pid() == 0) {
    if (i == 0) {
      fprintf (ferr, "i dt t ke\n");
      fp = fopen ("log", "w");
      fprintf(fp, "Capillary Scaling: Level %d tmax %g. Ohf %3.2e, Ec %3.2e, De infty\n", MAXlevel, tmax, Ohf, Ec);
      fprintf (fp, "i dt t ke\n");
      fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
      fclose(fp);
    } else {
      fp = fopen ("log", "a");
      fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
      fclose(fp);
    }
    fprintf (ferr, "%d %g %g %g\n", i, dt, t, ke);
  }

  assert(ke > -1e-10);
  // dump (file = "dump");
  if (ke < 1e-6 && i > 100){
    if (pid() == 0){
      fprintf(ferr, "kinetic energy too small now! Stopping!\n");
      fp = fopen ("log", "a");
      fprintf(fp, "kinetic energy too small now! Stopping!\n");
      fclose(fp);
    }
    return 1;
  }
}
