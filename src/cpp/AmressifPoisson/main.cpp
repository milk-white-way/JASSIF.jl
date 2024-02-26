// =================== LISTING KERNEL HEADERS ==============================
#include <AMReX_Gpu.H>
#include <AMReX_Utility.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_BCRec.H>
#include <AMReX_BCUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_VisMF.H>

#include "main.H"

// Modulization library
#include "fn_rhs.H"
#include "poisson.H"
#include "utilities.H"

using namespace amrex;

// ============================== MAIN SECTION ==============================//
/**
 * This is the code using AMReX for solving Navier-Stokes equation using
 * hybrid staggerred/non-staggered method
 * Note that the Contravariant variables stay at the face center
 * The pressure and Cartesian velocities are in the volume center
 */
int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    main_main();
    amrex::Finalize();
    return 0;
}

void main_main ()
{
    // What time is it now?  We'll use this to compute total run time.
    auto strt_time = ParallelDescriptor::second();

    // AMREX_SPACEDIM: number of dimensions
    // These are stock params for AMReX
    int n_cell, max_grid_size, nsteps, plot_int;

    // time = imported time from matlab
    Real time;

    // Physical boundary condition mapping
    // 0 is periodic
    // -1 is non-slip
    // 1 is slip
    Vector<int> bc_lo(AMREX_SPACEDIM, 0);
    Vector<int> bc_hi(AMREX_SPACEDIM, 0);

// =-=-=-=-=-=-=-=-=-=-=-=-=-=-= Parsing Inputs =-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    {
        // ParmParse is way of reading inputs from the inputs file
        ParmParse pp;

        // We need to get n_cell from the inputs file - this is the number of cells on each side of
        //   a square (or cubic) domain.
        pp.get("n_cell", n_cell);
        amrex::Print() << "INFO| number of cells in each side of the domain: " << n_cell << "\n";

        // The domain is broken into boxes of size max_grid_size
        pp.get("max_grid_size", max_grid_size);

        // Default plot_int to -1, allow us to set it to something else in the inputs file
        //  If plot_int < 0 then no plot files will be written
        plot_int = -1;
        pp.query("plot_int", plot_int);

        // Default nsteps to 1, allow us to set it to something else in the inputs file
        nsteps = 1;
        pp.query("nsteps", nsteps);

        // Reading current time in case of plot files
        pp.query("time", time);

        // read in BC
        pp.queryarr("bc_lo", bc_lo);
        pp.queryarr("bc_hi", bc_hi);
    }

    Vector<int> is_periodic(AMREX_SPACEDIM, 0);
    for (int idim=0; idim < AMREX_SPACEDIM; ++idim) {
        if (bc_lo[idim] == 0 && bc_hi[idim] == 0) {
            is_periodic[idim] = 1;
        }
        amrex::Print() << "INFO| periodicity in " << idim+1 << "th dimension: " << is_periodic[idim] << "\n";
    }

// =-=-=-=-=-=-=-=-=-=-=-=-=-=-= Defining System's Variables =-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // make BoxArray and Geometry
    BoxArray ba;
    Geometry geom;
    {
        IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
        IntVect dom_hi(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1));
        Box domain(dom_lo, dom_hi);

        // Initialize the boxarray "ba" from the single box "bx"
        ba.define(domain);
        // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
        ba.maxSize(max_grid_size);

        // Here, the real domain is a rectangular box defined by (0,0); (0,1); (1,0); (1,1)
        // This defines the physical box, [0,1] in each direction.
        RealBox real_box({AMREX_D_DECL( Real(0.0), Real(0.0), Real(0.0))},
                         {AMREX_D_DECL( Real(1.0), Real(1.0), Real(1.0))});

        // This defines a Geometry object
        //   NOTE: the coordinate system is Cartesian
        geom.define(domain, &real_box, CoordSys::cartesian, is_periodic.data());
    }

    // Nghost = number of ghost cells for each array
    int Nghost = 1; // 2nd order accuracy scheme is used for convective terms

    // Ncomp = number of components for userCtx
    int Ncomp = 2;

    // How Boxes are distrubuted among MPI processes
    // Distribution mapping between the processors
    DistributionMapping dm(ba);

    // UserContext MultiFab contains 2 components, pressure and Phi, at the cell center
    MultiFab userCtx(ba, dm, Ncomp, 0);
    MultiFab velCart(ba, dm, AMREX_SPACEDIM, 0);

    MultiFab poisson_rhs(ba, dm, 1, 0);
    MultiFab poisson_sol(ba, dm, 1, 0);

    //---------------------------------------------------------------
    // Defining the boundary conditions for each face of the domain
    // --------------------------------------------------------------
    Vector<BCRec> bc(poisson_sol.nComp());
    for (int n = 0; n < poisson_sol.nComp(); ++n)
    {
        for(int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            //Internal Dirichlet Periodic Boundary conditions, or bc_lo = bc_hi = 0
            if (bc_lo[idim] == BCType::int_dir) {
                bc[n].setLo(idim, BCType::int_dir);
            }
                //First Order Extrapolation for Neumann boundary conditions or bc_lo, bc_hi = 2
            else if (bc_lo[idim] == BCType::foextrap) {
                bc[n].setLo(idim, BCType::foextrap);
            }
                //External Dirichlet Boundary Condition, or bc_lo, bc_hi = 3
            else if(bc_lo[idim] == BCType::ext_dir) {
                bc[n].setLo(idim, BCType::ext_dir);
            }
            else {
                amrex::Abort("Invalid bc_lo");
            }

            //Internal Dirichlet Periodic Boundary conditions, or bc_lo = bc_hi = 0
            if (bc_hi[idim] == BCType::int_dir) {
                bc[n].setHi(idim, BCType::int_dir);
            }
                //First Order Extrapolation for Neumann boundary conditions or bc_lo, bc_hi = 2
            else if (bc_hi[idim] == BCType::foextrap) {
                bc[n].setHi(idim, BCType::foextrap);
            }
                //External Dirichlet Boundary Condition, or bc_lo, bc_hi = 3
            else if(bc_hi[idim] == BCType::ext_dir) {
                bc[n].setHi(idim, BCType::ext_dir);
            }
            else {
                amrex::Abort("Invalid bc_hi");
            }
        }
    }

    // Contravariant velocities live in the face center
    Array<MultiFab, AMREX_SPACEDIM> velCont;

    // gradient of phi
    Array<MultiFab, AMREX_SPACEDIM> grad_phi;

    // Due to the mismatch between the volume-center and face-center variables
    // The physical quantities living at the face center need to be blowed out one once in the respective direction
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
    {
        // flux(dir) has one component, zero ghost cells, and is nodal in direction dir
        BoxArray edge_ba = ba;
        edge_ba.surroundingNodes(dir);

        velCont[dir].define(edge_ba, dm, 1, 0);
        grad_phi[dir].define(edge_ba, dm, 1, 0);
    }

    // Print desired variables for debugging
    amrex::Print() << "INFO| number of dimensions: " << AMREX_SPACEDIM << "\n";
    amrex::Print() << "INFO| geometry: " << geom << "\n";
    amrex::Print() << "INFO| velCont[0] size: " << velCont[0].size() << "\n";

    amrex::Print() << "========================= IMPORT DATA STEP ========================= \n";
    // Current: Taylor-Green Vortex initial conditions
    // How partial periodic boundary conditions can be deployed?

    // Importing the data from the HDF5 file
    // The pressure field is temporarily stored in the poisson_sol MultiFab
    import_matlab_hdf5(poisson_sol, velCont, geom);

    //cont2cart(velCart, velCont, geom);

    //WriteSingleLevelPlotfile("pltPressure", poisson_sol, {"imp"}, geom, time, 0);

    // Transfer the pressure field to the userCtx MultiFab
    //MultiFab::Copy(userCtx, poisson_sol, 0, 0, 1, 0);
    // Reset the correction term - phi - to zero
    //poisson_sol.setVal(0.0);
    // Save phi to the userCtx MultiFab
    //MultiFab::Copy(userCtx, poisson_sol, 0, 1, 1, 0);
    
/*
    for (int n = 1; n <= nsteps; ++n)
    {
        amrex::Print() << "============================ ADVANCE STEP " << n << " ============================ \n";

        // Poisson solver
        //    Laplacian(\phi) = (Real(1.5)/dt)*Div(\hat{u}_i)

        // POISSON |1| Calculating the RSH
        poisson_righthand_side_calc(poisson_rhs, velCont, geom, dt);

        // POISSON |2| Init Phi at the begining of the Poisson solver
        poisson_sol.setVal(0.0);
        poisson_sol.FillBoundary(geom.periodicity());
        poisson_advance(poisson_sol, poisson_rhs, geom, ba, dm, bc);
        amrex::Print() << "SOLVING| finished solving Poisson equation. \n";
*/
        // Write a plotfile of the current data (plot_int was defined in the inputs file)
        //if (plot_int > 0)
        //{
        //    Export_Flow_Field("pltResults", userCtx, velCart, ba, dm, geom, time, 0);
        //}

        amrex::Print() << "========================== FINISH TIME: " << time << " ========================== \n";
/*
    }//end of time loop - this is the (n) loop!
*/

    // Call the timer again and compute the maximum difference
    // between the start time and stop time
    // over all processors
    auto stop_time = ParallelDescriptor::second() - strt_time;
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    ParallelDescriptor::ReduceRealMax(stop_time, IOProc);

    // Tell the I/O Processor to write out the "run time"
    amrex::Print() << "Run time = " << stop_time << std::endl;
}
