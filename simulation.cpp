/*
Below is a code that simulates the critical phenomena in the spherically symmetric'
gravitational collapse.

The code is written in C++ v11 using Qt v5.4 framework.

Please note that the usual naming styles/conventions for C++ coding, had been
sacrificed here, to match the mathematical notations and symbols that are used in
the supplementary PDF file.

Simulation raw data, are dumped to tab separated text files, in a directory that hard coded below.

To process, and generate various graphs of the generated raw data, see the Wolfram
Mathematica supplementary code.

Source repository, additional files, and licence can be found on:
https://github.com/TMS-Namespace/Critical-Phenomena-in-Gravitational-Collapse

*/

#define Includes {

    #include <QCoreApplication>
    #include <stdio.h>
    #include <stdlib.h>
    #include <math.h>
    #include <time.h>
    #include <iostream>
    #include <QFile>
    #include <QTextStream>
    #include <QStringBuilder>
    #include <QTime>
    #include <qmath.h>
    #include <limits.h>
    #include <sstream>

#define IncludesEnd }

//==========================================================

#define GlobalVars {

    #define GridSize 3000
    #define InvLn10 0.43429448190325182765112891891661
    #define Pi 3.1415926535897932384626433832795
    #define Separator "\t"

    // Represent a grid Point
    struct Point
    {
        int u;
        int v;
    };

    // Holds initial scalar field info
    struct InitialField
    {
        double Sigma;
        double v0;
        double A;
    };

    // This will hold a region info (-5,+5)
    struct Region
    {
        double Max;
        double Min;
    };

    double  a[GridSize][GridSize],
            Av[GridSize][GridSize],
            Phi[GridSize][GridSize],
            Fv[GridSize][GridSize],
            Fu[GridSize][GridSize],
            r[GridSize][GridSize],
            Ru[GridSize][GridSize],
            Rv[GridSize][GridSize],
            Error[GridSize][GridSize],
            T[GridSize];

    double h;
    /* Those are temp matrices, that will hold the value of the above matrices at
    some point*/
    double _Fv, _Fu, _Ru, _Rv, _Av, _r, _a;

#define GlobalVarsEnd }

//==========================================================

#define Prototypes {

    typedef double (*Solver)(Point, int, double);

    void    Evolve(int),
            Prepare_Variables(Point, int, double, double*),
            Prepare_Variables(Point, int, double),
            Integrate_Constraints(Point),
            Simulate(InitialField),
            Save_LnT_Lna(InitialField, double[GridSize], double), // "Ln" stands for logarithm
            Save_LnT_Phi(InitialField, double[GridSize], double),
            Save_Array(double[GridSize][GridSize], QString, int = 20),
            ProperTime(),
            Error_Calc(),
            MeasureTime(QString = ""),
            Save_LnT_LnR(InitialField, double[GridSize], double),
            Save_LnA_LnM(InitialField, Region, double);

    Region Find_Critical_Amplitude(InitialField, Region, double);

    bool IsVal(double);
    // IsCollapsing();

    QTextStream* Save_GetStream(QString);

    double  Solve_a(Point, int, double),
            Solve_r(Point, int, double),
            Solve_Rv(Point, int, double),
            Solve_Phi(Point, int, double),
            Solve_Ru(Point, int, double),
            Solve_Fu(Point, int, double),
            Solve_evl_Fv(double, double, double, double, double),
            Solve_evl_Av(double, double, double, double, double, double),
            Rung_Kutta_One_Cycle(Solver, Point),
            BH_Mass(Point), // BH stands for black hole
            Square(double),
            Ricci(Point);

    Point Collapsing_Point();

#define PrototypesEnd }

//==========================================================

int main(int argc, char* argv[])
{
	QCoreApplication ee(argc, argv);

	// QTime timer;
	// timer.start();

	// Increase console output precision to max
	std::cout.precision(std::numeric_limits<double>::digits10 + 1);

	// Scalar field parameters
	InitialField initial_field;
	initial_field.Sigma = 1 / qSqrt(2);
	initial_field.v0 = 2;

	//---------- Find Critical Amplitude
	Region region;
	region.Min = 0.05;
	region.Max = 0.06;
	Region A_Critical;
	// A_Critical = Find_Critical_Amplitude(initial_field, region, (region.Max - region.Min) / 10000);
	// Values for h = 0.0025, sigma = 1 / Sqrt(2), v0 = 2
	h = 0.0025;
	A_Critical.Max = 0.0513408202905885;
	A_Critical.Min = 0.0513408202905883;
	// Values for h = 0.0009, sigma = 1 / Sqrt(2), v0 = 2
	// h = 0.0009;
	// A_Critical.Max = 0.05134811985773;
	// A_Critical.Min = 0.05134811985766;

	// For v0 = 2, Sigma = 1 / qSqrt(2) BH forms approximately at v = 4.4982 thus minimum 
	// possible value for h = 4.4982 / 5000 = 0.0009 this Initial and Boundary conditions, 
	// critical amplitude is:
	// For h=0.0025 we get A * (0.051340, 0.051341)
	// For h=0.0020 we get A * (0.051342, 0.051343)
	// For h=0.0010 we get A * (0.051347, 0.051348)
	// For h=0.0009 we get A * (0.051348, 0.051349), A(0.05134811985766, 0.05134811985773)

	//----------Check if vacuum solution produces flat Minicowski space
	/*Simulate(initial_field);
	Save_Array(a, "a_Vacuum");
	Error_Calc();
	Save_Array(Error, "Error_Vacuum");*/

	//----------Save Subcritical results to text data files
	/*initial_field.A = A_Critical.Max - A_Critical.Max / 1000;
	Simulate(initial_field);

   Save_Array(a, "a_SuB");
   Save_Array(Phi, "Phi_SuB");
   Save_Array(r, "r_SuB");
   Save_Array(Rv, "Rv_SuB");
   Error_Calc();
   Save_Array(Error, "Error_SuB");*/

   //---------- Save Supercritical results
   /*initial_field.A = A_Critical.Max + 30 * A_Critical.Max / 1000; //*30 is just to make resultant graphics clearer
   Simulate(initial_field);
   Collapsing_Point();

   Save_Array(a, "a_SuP");
   Save_Array(Phi, "Phi_SuP");
   Save_Array(r, "r_SuP");
   Save_Array(Rv, "Rv_SuP");
   Error_Calc();
   Save_Array(Error, "Error_SuP");*/

   //============= Critical Behavior Validation

   //------1- We need to find ProperTime on Axis for Critical parameter
	initial_field.A = A_Critical.Max;
	Simulate(initial_field);
	Collapsing_Point(); // Just to Make sure that we have BH
	ProperTime();

	double T_Critical[GridSize];
	std::copy(std::begin(T), std::end(T), std::begin(T_Critical));

	//------2- Show self-similarity near critical parameter
	Save_LnT_Phi(initial_field, T_Critical, A_Critical.Max + A_Critical.Max / 100000000);
	Save_LnT_Phi(initial_field, T_Critical, A_Critical.Min - A_Critical.Min / 100000000);

	//------3- Show how a on axis changes with time
	/*Save_LnT_Lna(initial_field, T_Critical, A_Critical.Max + A_Critical.Max / 1000000);
	Save_LnT_Lna(initial_field, T_Critical, A_Critical.Min - A_Critical.Min / 1000000);*/

	//------ 4- Find how BH Mass changes
	initial_field.A = A_Critical.Max;
	// region.Max = A_Critical.Max + A_Critical.Max / 300000; region.Min = A_Critical.Max;
	// Save_LnA_LnM(initial_field, region, (region.Max - region.Min) / 30);

	//------ 5- Find Ricci Scalar
	/*Save_LnT_LnR(initial_field, T_Critical, A_Critical.Max + A_Critical.Max / 1000000);
	Save_LnT_LnR(initial_field, T_Critical, A_Critical.Min - A_Critical.Min / 1000000);*/

	std::cout << "Done..." << std::endl;
	return ee.exec();
}

//==========================================================

#define Core {

    inline Region Find_Critical_Amplitude(InitialField init_field, Region region, double accuracy)
    {

        // Those will hold info if BH formed on not
        bool from_BH, to_BH, middle_BH;

        std::cout << "Searching for Critical Amplitude in Interval A(" << region.Min << "," << region.Max << ")" << std::endl << std::endl;

        init_field.A = region.Min;
        Simulate(init_field);
        from_BH = (Collapsing_Point().u != -1);

        init_field.A = region.Max;
        Simulate(init_field);
        to_BH = (Collapsing_Point().u != -1);

        // Squeeze the region of critical amplitude from left and right
        while (region.Max - region.Min > accuracy) {
            std::cout << "Checking Interval A(" << region.Min << "," << region.Max << ")" << std::endl << " Fraction of target accuracy: " << (region.Max - region.Min) / accuracy << std::endl << std::endl;

            if (from_BH == to_BH) {
                region.Max = -1;
                region.Min = -1;
                return region;
            }
            else {
                init_field.A = (region.Max + region.Min) / 2;
                Simulate(init_field);
                middle_BH = (Collapsing_Point().u != -1);

                if (middle_BH == from_BH) {
                    region.Min = (region.Max + region.Min) / 2;
                }
                else {
                    region.Max = (region.Max + region.Min) / 2;
                    to_BH = middle_BH;
                }
            }
        }

        std::cout << "Critical Interval A(" << region.Min << "," << region.Max << ")" << std::endl;
        return region;
    }

    //==========================================================

    void ProperTime()
    {

        // Using Simpson rule
        int shift_from_axis{ 1 };

        T[0] = 0;
        T[1] = h * a[1][1 + shift_from_axis];

        for (int u = 2; u < GridSize - shift_from_axis; u++) {
            T[u] = T[u - 2] + h / 3 * (a[u - 1][u - 1 + shift_from_axis] + 4 * a[u][u + shift_from_axis] + a[u + 1][u + 1 + shift_from_axis]);
        }

        T[GridSize - 1] = T[GridSize - 2] + a[GridSize - 2][GridSize - 2 + shift_from_axis] * h;
    }

    //==========================================================

    // inline double BH_Mass(Point P){
    //    if (P.u > -1){
    //        return 0.5 * r[P.u][P.v] * (1 + 4 * Rv[P.u][P.v] * Ru[P.u][P.v]) / Sq(a[P.u][P.v]);
    //    }
    //    return 0;
    //}

    //==========================================================

    Point Collapsing_Point()
    {
        Point collaps_point;
        collaps_point.u = -1;
        collaps_point.v = -1;

        int MaxV{ -1 };

        for (int u = 0; u < GridSize; u++) {
            for (int v = u; v < GridSize; v++) {
                // Find where Rv changing it's sign
                if (Rv[u][v - 1] * Rv[u][v] < 0)
                {
                    std::cout << "Black Hole Found at v =" << v * h << std::endl;
                    // Within same slice, go backwards to find Max(r)
                    MaxV = v;

                    for (int BackV = v - 1; BackV > -1; BackV--) {
                        if (IsVal(r[u][BackV]) && r[u][BackV] > r[u][MaxV]) {
                            MaxV = BackV;
                        }
                    }

                    collaps_point.u = u;
                    collaps_point.v = MaxV;

                    std::cout << "r_Max at (u,v) = (" << u * h << "," << v * h << ")" << std::endl << std::endl;

                    return collaps_point;
                }
            }
        }

        // No Collapsing, Return (-1,-1)
        return collaps_point;
    }

    //==========================================================

    inline void Simulate(InitialField init_field)
    {

        std::cout << "Simulating GC for A=" << init_field.A << " v0 =" << init_field.v0 << " Sigma =" << init_field.Sigma << std::endl;

        // Set timer for calculation speed
        double calc_speed = 0;
        QTime timer;
        timer.start();

        // some Temporary values
        double diff, width, V;

        width = 1 / Square(init_field.Sigma);

        for (int v = 0; v < GridSize; v++) {
            V = v * h;
            Av[0][v] = 0;
            diff = V - init_field.v0;
            Fv[0][v] = 2 * init_field.A * V * (1 - V * diff * width) * exp(-width * Square(diff));
        }

        //--------------------------boundary conditions at the axis u=v

        a[0][0] = 1;                // Flat Space
        r[0][0] = h / 10;           // Minimum r to start from to avoid numerical problems
        Rv[0][0] = 0.5 * (a[0][0]); // Rv = a / 2 = -Ru
        Ru[0][0] = -Rv[0][0];
        Fu[0][0] = Fv[0][0];
        Phi[0][0] = 0;

        int u, v;

        Point P;
        P.u = 0;

        for (v = 0; v < GridSize - 1; v++) {
            P.v = v;
            Integrate_Constraints(P);
        }

        for (u = 1; u < GridSize; u++) {
            if (u % 100 == 0) {
                if (timer.elapsed() != 0)
                    calc_speed = 1000 * u / timer.elapsed();

                std::cout << "\33\r Slice: " << u << " Speed: " << calc_speed << " Slice/s";
            }

            // Boundary conditions
            for (v = 0; v < u; v++) {
                a[u][v] = 0;
                Phi[u][v] = 0;
                r[u][v] = 0;
                Rv[u][v] = 0;
                Ru[u][v] = 0;
                Av[u][v] = 0;
                Fu[u][v] = 0;
                Fv[u][v] = 0;
            }

            if (u == 1) {
                a[u][u] = a[u - 1][u + 1];
                Phi[u][u] = Phi[u - 1][u + 1];
            }
            else {
                a[u][u] = (4 * a[u - 1][u + 1] - a[u - 2][u + 2]) / 3;
                Phi[u][u] = (4 * Phi[u - 1][u + 1] - Phi[u - 2][u + 2]) / 3;
            }

            r[u][u] = h / 10; // Minimum r to start from to avoid numerical problems
            Rv[u][u] = 0.5 * (a[u][u]);
            Ru[u][u] = -Rv[u][u];
            Av[u][u] = Av[u - 1][u];
            Fv[u][u] = Fv[u - 1][u];
            Fu[u][u] = Fv[u][u];

            // Calculate the variables in the next u
            Evolve(u);
        }

        std::cout << "\33\r Slice: " << u << " Speed: " << calc_speed << " Slice/s" << std::endl << std::endl;
    }

    //==========================================================

    inline void Evolve(int u)
    {
        double a_Fv{ 0 }, a_Av{ 0 }, b_Ru{ 0 }, a_Rv{ 0 }, b_Rv{ 0 }, b_r{ 0 }, b_Fu{ 0 }, b_a{ 0 }, nominator{ 0 };

        Point P;

        for (int v = u + 1; v < GridSize - 1; v++) {

            // Evolution equations from u-1 --> u
            P.u = u - 1;
            P.v = v;
            Prepare_Variables(P, 1, 0);

            // E8
            a_Fv = _Fv + h * Solve_evl_Fv(_Ru, _Rv, _Fu, _Fv, _r);

            // M10
            a_Av = _Av + h * Solve_evl_Av(_Ru, _Rv, _Fu, _Fv, _r, _a);
            // Prepare_Variables(u, v-1, 1, 0);
            //============================

            // M14
            // Loading Solve_Rv will prepare the right variables for the next calculations too
            P.u = u;
            P.v = v - 1;
            a_Rv = _Rv + h * Solve_Rv(P, 1, 0);
            // Rv_ = _Rv + h * Solve_Rv(u, v-1, 1, 0);

            //============================ Explicit Solutions
            // B.1 solution of E2
            b_a = _a * (1 + 0.5 * h * _Av) / (1 - 0.5 * h * a_Av);

            // B.3 solution of E4
            b_r = _r + 0.5 * h * (_Rv + a_Rv);

            // B.4 solution of E14
            b_Rv = (_Rv + 0.5 * h * (2 * _Av * _Rv - _r * Square(_Fv) - b_r * Square(a_Fv))) / (1 - h * a_Av);

            nominator = 1 / (1 + h * b_Rv / (2 * b_r));

            // B.5 Solution of E11
            b_Ru = (_Ru - 0.5 * h * ((_Ru * _Rv + 0.25 * Square(_a)) / _r + Square(b_a) / (4 * b_r))) * nominator;

            // B.6 solution of E7
            b_Fu = (_Fu - 0.5 * h * ((_Ru * _Fv + _Rv * _Fu) / _r + b_Ru * a_Fv / b_r)) * nominator;

            // We need Fv/Av[u][v] in integration, thus we need to calc them
            //============================
            // Evolve q
            Fv[u][v] = 0.5 * (a_Fv + Fv[u - 1][v] + h * Solve_evl_Fv(b_Ru, b_Rv, b_Fu, a_Fv, b_r));

            // Evolve d
            Av[u][v] = 0.5 * (a_Av + Av[u - 1][v] + h * Solve_evl_Av(b_Ru, b_Rv, b_Fu, a_Fv, b_r, b_a));

            Integrate_Constraints(P);
        }
    }

    //==========================================================
    // Solve E8
    inline double Solve_evl_Fv(double Ru_, double Rv_, double Fu_, double Fv_, double r_)
    {
        return -(Ru_ * Fv_ + Rv_ * Fu_) / r_;
    }

    //==========================================================
    // Solve E10
    inline double Solve_evl_Av(double Ru_, double Rv_, double Fu_, double Fv_, double r_, double a_)
    {
        return (Ru_ * Rv_ + 0.25 * Square(a_)) / Square(r_) - Fv_ * Fu_;
    }

    //==========================================================

    inline void Integrate_Constraints(Point P)
    {
        int V = P.v + 1, u = P.u, v = P.v;

        a[u][V] = a[u][v] + Rung_Kutta_One_Cycle(Solve_a, P);

        r[u][V] = r[u][v] + Rung_Kutta_One_Cycle(Solve_r, P);

        Rv[u][V] = Rv[u][v] + Rung_Kutta_One_Cycle(Solve_Rv, P);

        Phi[u][V] = Phi[u][v] + Rung_Kutta_One_Cycle(Solve_Phi, P);

        Ru[u][V] = Ru[u][v] + Rung_Kutta_One_Cycle(Solve_Ru, P);

        Fu[u][V] = Fu[u][v] + Rung_Kutta_One_Cycle(Solve_Fu, P);
    }

#define CoreEnd}

//==========================================================

#define	Helpers{

    inline void Error_Calc()
    {
        MeasureTime();

        std::cout << "Calculating Error Matrix..." << std::endl;

        Point P;

        // Checking accuracy of equation E12
        for (int u = 2; u < GridSize - 2; u++) {
            for (int v = u; v < GridSize - 2; v++) {

                // using 8 order accuracy centric difference method to calc D[Rv,u]
                double Rv_u;
                Rv_u = (45 * (Rv[u + 1][v] - Rv[u - 1][v]) - 9 * (Rv[u + 2][v] - Rv[u - 2][v]) + (Rv[u + 3][v] - Rv[u - 3][v])) / (60 * h);

                P.u = u;
                P.v = v;
                Prepare_Variables(P, 1, 0);

                Error[u][v] = _r * Rv_u + _Rv * _Ru + 0.25 * Square(_a);
            }
        }

        MeasureTime("Error Matrix Calculation");
    }

    //==========================================================

    inline double Square(double v)
    {
        return v * v;
    }

    //==========================================================

    QTime Timer;

    inline void MeasureTime(QString operation_name)
    {
        if (operation_name == "") {
            Timer.start();
        }
        else {
            std::cout << operation_name.toLocal8Bit().constData() << " time: " << Timer.elapsed() / (1000) << " sec." << std::endl
                << std::endl;
        }
    }

    //==========================================================

    inline double Ricci(Point P)
    {
        int u = P.u, v = P.v;
        return -32 * Pi * Fu[u][v] * Fv[u][v] / Square(a[u][v]);
    }

    //==========================================================

    /* This will check if the value is a reasonable number and not NAN or Inf, so that we
    will not print it out, because Wolfram Mathemtica can't deal with it. Warning: this uses lot
    of CPU recourses.*/
    inline bool IsVal(double val)
    {
        QString str_val = QString::number(val);

        if (val != val || str_val == "nan") {
            // return "NaN";
            return false;
        }
        else if (val > std::numeric_limits<qreal>::max()) {
            // return "+Inf";
            return false;
        }
        else if (val < -std::numeric_limits<qreal>::max()) {
            // return "-Inf";
            return false;
        }
        else {
            return true;
        }
    }

#define HelpersEnd}

//==========================================================

#define DumpingData{

    void Save_LnT_LnR(InitialField init_critical_field, double T_Critical[GridSize], double A_For)
    {
        MeasureTime();

        QTextStream* file_stream;

        if (init_critical_field.A < A_For) {
            file_stream = Save_GetStream("LnT_LnR_Sup");
        }
        else {
            file_stream = Save_GetStream("LnT_LnR_Sub");
        }

        init_critical_field.A = A_For;
        Simulate(init_critical_field);
        ProperTime();

        double t, f;
        Point P;

        for (int u = 0; u < GridSize; u++) {
            t = qLn(qAbs(T_Critical[u] - T[u]));
            P.u = u;
            P.v = u;
            f = qLn(1 - Ricci(P));

            if (IsVal(t) && IsVal(f)) {
                (*file_stream) << t << Separator << f << endl;
            }
        }

        file_stream->flush();

        MeasureTime("Ln(T*-T)-->Ln(1-R)");
        (*file_stream).~QTextStream();
    }

    //==========================================================

    void Save_LnA_LnM(InitialField init_critical_field, Region region, double A_Step)
    {
        MeasureTime();

        QTextStream* file_stream = Save_GetStream("LnA_LnM");

        int step{ 1 };
        double A{ 0 }, A_Critical = init_critical_field.A;

        for (A = region.Min; A <= region.Max; A = A + A_Step) {
            std::cout << "--------------- Step " << step << " of " << (int)((region.Max - region.Min) / A_Step) << " Steps." << std::endl << std::endl;
            step++;

            init_critical_field.A = A;
            Simulate(init_critical_field);
            Point P = Collapsing_Point();

            if (P.v == -1 || P.u == -1) {
                std::cout << "No BH formation point found for A = " << A << "." << std::endl;
                return;
            }

            (*file_stream) << P.u << Separator << P.v << Separator << qLn(A - A_Critical) << Separator << qLn(r[P.u][P.v] / 2) << endl;
        }

        MeasureTime("Ln(A)-->Ln(M)");
        (*file_stream).~QTextStream();
    }

    //==========================================================

    void Save_LnT_Lna(InitialField init_critical_field, double T_Critical[GridSize], double A_For)
    {
        MeasureTime();
        QTextStream* file_stream;

        if (A_For > init_critical_field.A) {
            file_stream = Save_GetStream("LnT_Lna_SuP");
        }
        else {
            file_stream = Save_GetStream("LnT_Lna_SuB");
        }

        init_critical_field.A = A_For;
        Simulate(init_critical_field);
        Collapsing_Point(); // Just to Make sure that we have BH
        ProperTime();

        double lt{ 0 }, la{ 0 };

        for (int u = 0; u < GridSize; u++) {
            lt = qLn(qAbs(T_Critical[u] - T[u]));
            la = qLn(a[u][u + 1]);

            if (IsVal(lt) && IsVal(la))
                (*file_stream) << lt << Separator << la << endl;
        }

        MeasureTime("Ln(T*-T)-->Ln(Phi)");
        (*file_stream).~QTextStream();
    }

    //==========================================================

    void Save_LnT_Phi(InitialField init_critical_field, double T_Critical[GridSize], double A_For)
    {
        MeasureTime();

        QTextStream* file_stream;

        if (init_critical_field.A < A_For) {
            file_stream = Save_GetStream("LnT_Phi_Sup");
        }
        else {
            file_stream = Save_GetStream("LnT_Phi_Sub");
        }

        init_critical_field.A = A_For;
        Simulate(init_critical_field);
        ProperTime();

        double t, f;

        for (int u = 0; u < GridSize; u++) {
            t = qLn(qAbs(T_Critical[u] - T[u]));
            f = Phi[u][u];

            if (IsVal(t) && IsVal(f)) {
                (*file_stream) << t << Separator << f << endl;
            }
        }

        file_stream->flush();

        MeasureTime("Ln(T*-T)-->Phi");
        (*file_stream).~QTextStream();
    }

    //==========================================================

    QTextStream* Save_GetStream(QString file_name)
    {
        QString file_path = "C:\\Simulation Results\\" + file_name + ".cvs";

        if (QFile::exists(file_path))
            QFile::remove(file_path);

        QFile* file = new QFile(file_path);

        if ((*file).open(QIODevice::ReadWrite)) {
            std::cout << "Dumping Results to " << file_name.toLocal8Bit().constData() << ".cvs" << std::endl;

            // QTextStream s(&File);
            // SaveStream = new QTextStream(&File);
            QTextStream* file_stream = new QTextStream(File);
            // (*file_stream).setRealNumberNotation(QTextStream::FixedNotation);
            // (*file_stream).setRealNumberPrecision(20);
            return file_stream;
        }
    }

    //==========================================================

    inline void Save_Array(double array[GridSize][GridSize], QString file_name, int point_step)
    {
        MeasureTime();
        QTextStream* file_stream = Save_GetStream(file_name);

        double val;

        for (int u = 0; u < GridSize; u = u + point_step) {
            for (int v = u; v < GridSize; v = v + point_step) {
                val = array[u][v];

                if (IsVal(val) && qAbs(val) < 5)
                {
                    (*file_stream) << u << Separator << v << Separator << val << endl;
                }
                else
                {
                    (*file_stream) << u << Separator << v << Separator << -1 << endl;
                }
            }
        }

        MeasureTime("[u,v]-->" + file_name);
        (*file_stream).~QTextStream();
    }

#define DumpingDataEnd}

//==========================================================

#define Solvers {

    // Solving E2
    inline double Solve_a(Point P, int K_level, double previous_K)
    {
        Prepare_Variables(P, K_level, previous_K, &_a);
        return _a * _Av;
    }

    // Solving E4
    inline double Solve_r(Point P, int K_level, double previous_K)
    {
        /*The problem in calculating r, is that we need Rv[u][v+1] when calculating k2 in
        Prepare_Variables() function, which supposed to be calculated by Rung-Kutta method
        only latter, to avoid that, and as an exception, we will calculate it here explicitly
        from equation E14 by Euler method, anyway, it will be calculated and overwritten by
        a more precise value by Rung-Kutta in the next solver..
        Also, we need to calculate it just once for same u,v..*/
        int u = P.u, v = P.v;

        if (K_level == 1) {
            Rv[u][v + 1] = Rv[u][v] + h * Solve_Rv(P, 1, 0);
        }
        else {
            // now calling Prepare_Variables will give the right result..
            Prepare_Variables(P, K_level, previous_K);
        }

        return _Rv;
    }

    // Solving E14
    inline double Solve_Rv(Point P, int K_level, double previous_K)
    {
        // Because Y here is _Rv, it will not be affected by calculation done in Solve_r
        Prepare_Variables(P, K_level, previous_K, &_Rv);
        return 2 * _Av * _Rv - _r * Square(_Fv);
    }

    // Solving E6
    inline double Solve_Phi(Point P, int K_level, double previous_K)
    {
        Prepare_Variables(P, K_level, previous_K);
        return _Fv;
    }

    // Solving E11
    inline double Solve_Ru(Point P, int K_level, double previous_K)
    {
        Prepare_Variables(P, K_level, previous_K, &_Ru);
        return -(_Ru * _Rv + Square(_a) * 0.25) / _r;
    }

    // Solving E7
    inline double Solve_Fu(Point P, int K_level, double previous_K)
    {
        Prepare_Variables(P, K_level, previous_K, &_Fu);
        return -(_Ru * _Fv + _Rv * _Fu) / _r;
    }

#define SolversEnd }

//==========================================================

#define RungKuttaHelpers {

    /*For Some solvers, there is no Y at all in G of the equation dY/dv=G(Y,v)=G(v), but
    because there is no way to assign a default value for optional parameter of pointer
    type in a function in C++, we need to overload/wrap it.*/
    inline void Prepare_Variables(Point P, int K_level, double previous_K)
    {
        Prepare_Variables(P, K_level, previous_K, NULL);
    }

    //==========================================================

    inline void Prepare_Variables(Point P, int K_level, double previous_K, double* Y)
    {
        int V = P.v, u = P.u, v = P.v;

        if (K_level > 1)
            V = v + 1;

        // Phi, Au are not used directly in any of the equations
        _Av = Av[u][V];

        _Ru = Ru[u][V];
        _Rv = Rv[u][V];

        _Fv = Fv[u][V];
        _Fu = Fu[u][V];

        _r = r[u][V];
        _a = a[u][V];

        if (K_level == 2 || K_level == 3) {
            _Av = (Av[u][v] + _Av) * 0.5;

            _Ru = (Ru[u][v] + _Ru) * 0.5;
            _Rv = (Rv[u][v] + _Rv) * 0.5;

            _Fv = (Fv[u][v] + _Fv) * 0.5;
            _Fu = (Fu[u][v] + _Fu) * 0.5;

            _r = (r[u][v] + _r) * 0.5;
            _a = (a[u][v] + _a) * 0.5;
        }

        if (Y != NULL && K_level > 1)
        {
            /* In case k4 we add k3 instead of k3 / 2 as in others, so we multiply it by 2
            to unify equations*/
            if (K_level == 4)
                previous_K = previous_K * 2;

            if (Y == &_Fu)
                _Fu = Fu[u][v] + previous_K * 0.5;
            if (Y == &_Ru)
                _Ru = Ru[u][v] + previous_K * 0.5;
            if (Y == &_Rv)
                _Rv = Rv[u][v] + previous_K * 0.5;
            if (Y == &_a)
                _a = a[u][v] + previous_K * 0.5;
        }
    }

    //==========================================================

    inline double Rung_Kutta_One_Cycle(Solver f, Point P)
    {
        double k1{ 0 }, k2{ 0 }, k3{ 0 }, k4{ 0 };

        k1 = h * f(P, 1, 0);
        k2 = h * f(P, 2, k1);
        k3 = h * f(P, 3, k2);
        k4 = h * f(P, 4, k3);

        return (k1 + k4 + 2 * (k2 + k3)) / 6;
    }

#define RungKuttaHelpersEnd }