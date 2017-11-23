#include <acado_toolkit.hpp>
#include <acado_gnuplot/gnuplot_window.hpp>

int main( ) {
    USING_NAMESPACE_ACADO
            
    DifferentialState x1;
    DifferentialState x2;
    DifferentialState x3;
    DifferentialState x4;
    DifferentialState x5;
    
    Control U;
    
    double a11= 0.9988;
    double a12= -0.06006;
    double a13= 0.03115;
    double a14= -0.001041;
    double a15= -0.00007769;
    double a21= 0.06216;
    double a22= 0.9947;
    double a23= -0.07866;
    double a24= 0.002813;
    double a25= 0.0002063;    
    double a31= -0.02486;
    double a32= 0.0794;
    double a33= 0.9928;
    double a34= -0.0844;
    double a35= -0.004281;
    double a41= -0.004255;
    double a42= 0.002692;
    double a43= 0.0803;
    double a44= 0.9875;
    double a45= 0.1455;
    double a51= -0.001548;
    double a52= -0.0002731;
    double a53= -0.007124;
    double a54= -0.1251;
    double a55= 0.8611;
    
    double b11= -0.01064;
    double b21= -0.05861;
    double b31= -0.2075;
    double b41= -0.3704;
    double b51= 0.1243;
    
    double c11= -1.731;
    double c12= -0.5872;
    double c13= -0.003244;

    // Defina una ecuación diferencial
    DifferentialEquation f;
    
    f << dot(x1) == a11*x1 + a12*x2 + a13*x3 + a14*x4 + a15*x5 + b11*U;
    f << dot(x2) == a21*x1 + a22*x2 + a23*x3 + a24*x4 + a25*x5 + b21*U;
    f << dot(x3) == a31*x1 + a32*x2 + a33*x3 + a34*x4 + a35*x5 + b31*U;
    f << dot(x4) == a41*x1 + a42*x2 + a43*x3 + a44*x4 + a45*x5 + b41*U;
    f << dot(x5) == a51*x1 + a52*x2 + a53*x3 + a54*x4 + a55*x5 + b51*U;
    
    // Defina Función minimos cuadrados
    Function h;
    
    h << x1;
    h << x2;
    h << x3;
    h << x4;
    h << x5;
    h << U;

    Function g;
    
    g << c11*x1 + c12*x2 + c32*x3;
    
    Matrix Q(6,6);
    Q(0,0) = 10.0;
    Q(1,1) = 10.0;
    Q(2,2) = 10.0;
    Q(3,3) = 1.0;
    Q(4,4) = 1.0;
    Q(5,5) = 1.0;
    
    // Referencia
     Vector r(0);
     r.setAll( 27.0 );
    
    // Defina objeto de problema de control optimo
    const double tStart = 0.0;
    const double tEnd = 3.0;
    
    OCP.ocp(tStart, tEnd, 20);
    
    ocp.minimizeLSQ( Q, g, r );
    
    ocp.subjectTo( f );
    
    ocp.subjectTo(0.0 <= U <= 6.0 );
    
    // Algoritmo tiempo Real
    
    RealTimeAlgorithm alg ( ocp, 0.025 );
    alg.set( MAX_NUM_ITERATIONS, 1 );
    alg.set( PLOT_RESOLUTION, MEDIUM );
    
    GnuplotWindow window;
        window.addSubplot(x1, "State I" );
        window.addSubplot(x2, "State II" );
        window.addSubplot(x3, "State III" );
        window.addSubplot(x4, "State IV" );
        window.addSubplot(x5, "State V" );
        window.addSubplot(U, "Control" );
        window.addSubplot(g, "Salida" );
        alg << window;
    
    // Controlador  REFERENCIA TCP IP
        
   StaticReferenceTrajectory zeroReference( 27.0 );
    
   Controller controller( alg, zeroReference );
    
     Vector y( 0 );
     y(0) = 20.0;
    
    controller.init( 0.0,y );
    controller.step( 0.0,y );
    
    // https://www.youtube.com/watch?v=is11f5ac6lA
	// MHE PROBLEM FORMULATION
	//
	OCP ocp1(0.0, N * Ts, N);

	ocp1.subjectTo( f );

	// Measurement function h(x, u) on first N nodes
 	Function hm; // Medida TCP  IP
 
 	hm << x1 << x2 << x3 << x4 << x5 << U;

	// Weighting matrices and measurement functions
	DMatrix W = eye<double>( 6 );


	Function hN;
	hN << x1 << x2 << x3 << x4 << x5 << U;

	DMatrix WN = eye<double>( 6 );
	WN(0, 0) = W(0, 0);
	WN(1, 1) = W(1, 1);
	WN(2, 2) = W(2, 2);
	WN(3, 3) = W(3, 3);
	WN(4, 4) = W(4, 4);
    WN(5, 5) = W(5, 5);
	ocp1.minimizeLSQ(W, hm);
	ocp1.minimizeLSQEndTerm(WN, hN);

	OCPexport mhe( ocp1 );

	mhe.set(INTEGRATOR_TYPE, INT_RK4);
	mhe.set(NUM_INTEGRATOR_STEPS, N * Ni);

	mhe.set(HESSIAN_APPROXIMATION, GAUSS_NEWTON);
	mhe.set(DISCRETIZATION_TYPE, MULTIPLE_SHOOTING);

	mhe.set(HOTSTART_QP, YES);
// NOTE: This is crucial for export of MHE!
	mhe.set(SPARSE_QP_SOLUTION, CONDENSING);
	mhe.set(FIX_INITIAL_STATE, NO);

//	mhe.set( LEVENBERG_MARQUARDT, 1e-10 );

	if (mhe.exportCode("Lsd_mhe_export") != SUCCESSFUL_RETURN)
		exit( EXIT_FAILURE );

	mhe.printDimensionsQP( );
     
    return 0;
    
}