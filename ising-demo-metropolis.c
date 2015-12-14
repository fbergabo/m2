#include <stdio.h>
#include <math.h>

#include "ising_datastructure.h"
#include "ising.h"
#include "metropolis.h"
#include "progressbar.h"
#include "ising-demo-metropolis.h"
#include "plot_fun.c"
#define STEPS 20                //globally defined number of steps in termperature loop

int main (int argc, char **argv)
{

    int nt = STEPS;

    //arrays to store temperature and specific heat for each step
    double spec[STEPS], temp[STEPS];

    struct ising model;

    //define starting and ending temperature, and step size
    double tmin = 3., tmax = 4., dt, t = tmin;
    double h = 0.;              // magnetic field
    double j = 1.;              // ferromagnetic coupling

    dt = (tmax - tmin) / (nt - 1);
    int lx = 32;
    int ly = lx;

    //saves temperatures into temp[] array
    for (int bg = 0; bg < nt; bg++)
    {
        temp[bg] = tmin + bg * dt;
    }

    int mcsteps = 1024 * 128;   // 2^{17}

    //contingency for inability to allocate spins
    if (ising_init (&model, lx, ly, j, h, t) != 0)
    {
        fprintf (stderr, "%s - cannot allocate spin array\n", argv[0]);
        argc = 1;
        exit (1);
    }

    //prints the header
    print_header (&model, stderr);

    int therm_steps = mcsteps / 4;
    int report = therm_steps / 32;

    //loops for number of steps 
    for (int q = 0; q < nt; q++)
    {

        //loops over thermalization steps
        for (int s = 1; s <= therm_steps; s++)
        {
            //thermalization
            (void) one_metropolis_step_per_spin (&model);
            if (s % report == 0)
            {
                text_progress_bar (s, therm_steps, 32, "Thermalization",
                    stderr);
            }
        }
        fprintf (stderr, " Done\n");


        report = mcsteps / 32;

        //initialize placeholders
        double mav = 0., m2av = 0., eav = 0., e2av = 0.;
        double cv, xi;


        for (int s = 1; s <= mcsteps; s++)
        {
            //production
            (void) one_metropolis_step_per_spin (&model);
            double m = magnetization_per_spin (&model);
            double e = energy_per_spin (&model);

            mav += m;
            m2av += m * m;
            eav += e;
            e2av += e * e;

            if (s % report == 0)
            {
                text_progress_bar (s, mcsteps, 32, "Production", stderr);
            }

        }
        fprintf (stderr, " Done\n\n");

        mav /= mcsteps;
        eav /= mcsteps;
        m2av /= mcsteps;
        e2av /= mcsteps;

        // specific heat dE/dt
        cv = (e2av - eav * eav) / (t * t);
        spec[q] = cv;
        // magnetic susceptibility
        xi = (m2av - mav * mav) / t;
        printf ("%.3f   %.3f  % f %f  % f   %f\n", t, h, mav, xi, eav, cv);
        t += dt;
        ising_reinit (&model, h, t);
    }
    ising_free (model);
    plot (temp, spec, nt);

    return 0;
}

//Header function
void print_header (struct ising *model, FILE * fp)
{
    fprintf (fp, " Two-dimensional Ising model - Metropolis simulation\n");
    fprintf (fp, " ---------------------------------------------------\n");
    fprintf (fp, " T/J = %.3f  H/J = %.3f  %dx%d\n\n", model->T / model->J,
        model->H / model->J, model->Lx, model->Ly);
}
