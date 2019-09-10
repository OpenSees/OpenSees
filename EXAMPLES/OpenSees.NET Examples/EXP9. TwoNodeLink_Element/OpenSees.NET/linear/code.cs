public static void Exp9_linear()
        {
            var theDomain = new OpenSees.Components.DomainWrapper();
            var Tension = 10000.0 / 5.0;
            var IP = 10000.0 / 5.0;
            var OOP = 10000.0 / 5.0;
            var about_Tension = 10000.0 / 5.0;
            var about_IP = 10000.0 / 5.0;
            var about_OOP = 10000.0 / 5.0;

            var E_Timber_Long = 13200.0;
            var G_Timber_Long = 240.0;

            var tickness = 40.0;

            // materials
            var Stiffness_Tension = 10;
            var Stiffness_IP = 11;
            var Stiffness_OOP = 12;
            var Stiffness_about_Tension = 13;
            var Stiffness_about_IP = 14;
            var Stiffness_about_OOP = 15;

            var Stiffness_TensionMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(Stiffness_Tension, Tension, 0);
            //uniaxialMaterial Elastic	$Stiffness_Tension[expr	$Tension]

            var Stiffness_IPMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(Stiffness_IP, IP, 0);
            //uniaxialMaterial Elastic	$Stiffness_IP[expr	$IP]

            var Stiffness_OOPMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(Stiffness_OOP, OOP, 0);
            //uniaxialMaterial Elastic	$Stiffness_OOP[expr	$OOP]

            var Stiffness_about_TensionMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(Stiffness_about_Tension, about_Tension, 0);
            //uniaxialMaterial Elastic	$Stiffness_about_Tension[expr	$about_Tension]

            var Stiffness_about_IPMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(Stiffness_about_IP, about_IP, 0);
            //uniaxialMaterial Elastic	$Stiffness_about_IP[expr	$about_IP]

            var Stiffness_about_OOPMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(Stiffness_about_OOP, about_OOP, 0);
            //uniaxialMaterial Elastic	$Stiffness_about_OOP[expr	$about_OOP]

            
            /*
            node	1			0.0			0.0			0.0
            node	2			100.0		-250.0		600.0
             */

            theDomain.AddNode(new OpenSees.Components.NodeWrapper[] {
                new OpenSees.Components.NodeWrapper(1,6,0,0,0),
                new OpenSees.Components.NodeWrapper(2,6,10,-10,40),
            });

            //fix	1		1	1	1	1	1	1
            for (var i = 0; i < 6; i++)
                theDomain.AddSP_Constraint(
                    new OpenSees.Components.Constraints.SP_ConstraintWrapper(1, i, 0, true)
                    );

            //set Node_1_x            0.
            //set Node_1_y            0.
            //set Node_1_z            0.

            //set Node_2_x            100.0
            //set Node_2_y            -250.0
            //set Node_2_z            600.0


            //set Vec_X_x[expr   $Node_2_x -$Node_1_x]
            //set Vec_X_y[expr   $Node_2_y -$Node_1_y]
            //set Vec_X_z[expr   $Node_2_z -$Node_1_z]

            //set Vec_Y_x[expr -$Vec_X_y]
            //set Vec_Y_y[expr   $Vec_X_x]
            //set Vec_Y_z             0.0

            //set Vec_Z_x[expr   $Vec_X_y *$Vec_Y_z -$Vec_X_z *$Vec_Y_y]
            //set Vec_Z_y[expr   $Vec_X_z *$Vec_Y_x -$Vec_X_x *$Vec_Y_z]
            //set Vec_Z_z[expr   $Vec_X_x *$Vec_Y_y -$Vec_X_y *$Vec_Y_x]


            //puts    $Vec_Z_x
            //puts	$Vec_Z_y
            //puts	$Vec_Z_z



            //set Trans_tag   1
            //geomTransf Linear $Trans_tag    $Vec_Z_x    $Vec_Z_y    $Vec_Z_z

            var Node_1_x = 0.0;
            var Node_1_y = 0.0;
            var Node_1_z = 0.0;

            var Node_2_x = 100.0;
            var Node_2_y = -250.0;
            var Node_2_z = 600.0;


            var Vec_X_x = Node_2_x - Node_1_x;
            var Vec_X_y = Node_2_y - Node_1_y;
            var Vec_X_z = Node_2_z - Node_1_z;

            var Vec_Y_x = -Vec_X_y;
            var Vec_Y_y = Vec_X_x;
            var Vec_Y_z = 0.0;

            var Vec_Z_x = Vec_X_y * Vec_Y_z - Vec_X_z * Vec_Y_y;
            var Vec_Z_y = Vec_X_z * Vec_Y_x - Vec_X_x * Vec_Y_z;
            var Vec_Z_z = Vec_X_x * Vec_Y_y - Vec_X_y * Vec_Y_x;

            
            //element twoNodeLink		1		1		2		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	$Vec_Y_x	$Vec_Y_y	$Vec_Y_z
            var ele = new OpenSees.Elements.TwoNodeLinkWrapper(1, 3, 1, 2, new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 }),
                new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat },
                new OpenSees.VectorWrapper(new double[] { Vec_Y_x, Vec_Y_y, Vec_Y_z }), new OpenSees.VectorWrapper(0),
                new OpenSees.VectorWrapper(/*new double[] { 0.5, 0.5, 0.5, 0.5 }*/0), new OpenSees.VectorWrapper(0), 0, 0);

            theDomain.AddElement(ele);


            //            pattern Plain 1 Linear {

            //                load    2       100.        200.        300.        50.     100.        50.

            //}
            var theSeries = new OpenSees.Components.Timeseries.LinearSeriesWrapper();
            var theLoadPattern = new OpenSees.Components.LoadPatterns.LoadPatternWrapper(1);
            theLoadPattern.SetTimeSeries(theSeries);
            theDomain.AddLoadPattern(theLoadPattern);

            var theLoadValues = new OpenSees.VectorWrapper(new double[] { 100, 200, 300, 50, 100, 50 });
            theDomain.AddNodalLoad(new OpenSees.Components.Loads.NodalLoadWrapper(1, 2, theLoadValues, false), 1);

            //recorder Node -file	$folder/Displacements/Node_2_Disp.out	-node 2 	-dof	1	2	3	4	5	6	disp

            var savepath = System.Environment.CurrentDirectory +
                @"\100_Results_exp9_linear\Displacements\";
            if (!System.IO.Directory.Exists(savepath))
                System.IO.Directory.CreateDirectory(savepath);

            {
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var nodeTags = new int[] { 2 };
                foreach (var tag in nodeTags)
                {
                    var opsstream = new OpenSees.Handlers.DataFileStreamWrapper(savepath + $@"\Node_{tag}_Disp.out");
                    var recorder = new OpenSees.Recorders.NodeRecorderWrapper(dirIds, new OpenSees.IDWrapper(new int[] { tag }), 0, "disp", theDomain, opsstream);
                    theDomain.AddRecorder(recorder);
                }
            }


            // analyze

            var Dincr = 0.01;
            var Dmax = 1.0 / Dincr;

            var theModel = new OpenSees.AnalysisModelWrapper();

            var theSolnAlgo = new OpenSees.Algorithms.KrylovNewtonWrapper();

            var theIntegrator = new OpenSees.Integrators.Static.LoadControlWrapper
                (Dincr, 1, Dincr, Dincr);

            var theHandler = new OpenSees.Handlers.TransformationConstraintHandlerWrapper();

            var theRCM = new OpenSees.GraphNumberers.RCMWrapper(false);
            var theNumberer = new OpenSees.Numberers.DOF_NumbererWrapper(theRCM);

            var theSolver = new OpenSees.Systems.Linears.BandGenLinLapackSolverWrapper();
            var theSOE = new OpenSees.Systems.Linears.BandGenLinSOEWrapper(theSolver);

            var Tol = 1.0e-6;                        // Convergence Test: tolerance
            var maxNumIter = 400;                // Convergence Test: maximum number of iterations that will be performed before "failure to converge" is returned
            var printFlag = 0;                // Convergence Test: flag used to print information on convergence (optional)        # 1: print information on each step; 
                                              //var TestType EnergyIncr;	# Convergence-test type

            var test = new OpenSees.ConvergenceTests.CTestEnergyIncrWrapper
                (Tol, maxNumIter, printFlag, 2);
            var theAnalysis = new OpenSees.Analysis.StaticAnalysisWrapper(
                theDomain,
                theHandler,
                theNumberer,
                theModel,
                theSolnAlgo,
                theSOE,
                theIntegrator,
                test);

            theAnalysis.Analyze((int)Dmax);

            Console.WriteLine("done....");
            Console.ReadKey();

        }