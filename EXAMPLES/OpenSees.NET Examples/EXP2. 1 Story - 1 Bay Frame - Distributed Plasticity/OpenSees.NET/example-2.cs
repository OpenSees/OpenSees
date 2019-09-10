		public static void Exp2()
        {
            var theDomain = new OpenSees.Components.DomainWrapper();
            var e = 2.0e11;
            var i = 6.572e-5;
            var a = 4.265e-3;

            var nodes = new OpenSees.Components.NodeWrapper[]
            {
                new OpenSees.Components.NodeWrapper(1,3,0,0),
                new OpenSees.Components.NodeWrapper(2,3,0,3),
                new OpenSees.Components.NodeWrapper(3,3,0,3),
                new OpenSees.Components.NodeWrapper(4,3,4,3),
                new OpenSees.Components.NodeWrapper(5,3,4,0),
                new OpenSees.Components.NodeWrapper(6,3,2,3),
            };

            theDomain.AddNode(nodes);

            var sps = new OpenSees.Components.Constraints.SP_ConstraintWrapper[]
            {
                new OpenSees.Components.Constraints.SP_ConstraintWrapper(1, 0, 0.0, true),
                new OpenSees.Components.Constraints.SP_ConstraintWrapper(1, 1, 0.0, true),
                new OpenSees.Components.Constraints.SP_ConstraintWrapper(1, 2, 0.0, true),
                new OpenSees.Components.Constraints.SP_ConstraintWrapper(5, 0, 0.0, true),
                new OpenSees.Components.Constraints.SP_ConstraintWrapper(5, 1, 0.0, true),
                new OpenSees.Components.Constraints.SP_ConstraintWrapper(5, 2, 0.0, true),
            };
            theDomain.AddSP_Constraint(sps);


            //equalDOF 2 3 1 2
            var Ccr = new OpenSees.MatrixWrapper(2, 2);
            Ccr.Zero();
            Ccr[0, 0] = 1.0;
            Ccr[1, 1] = 1.0;
            var rcDOF = new OpenSees.IDWrapper(2);
            rcDOF[0] = 0;
            rcDOF[1] = 1;
            theDomain.AddMP_Constraint(new OpenSees.Components.Constraints.MP_ConstraintWrapper(2, 3, Ccr, rcDOF, rcDOF));


            var transftagCrdT = new OpenSees.Elements.CrdTransfs.LinearCrdTransf2dWrapper();
            var ele1 = new OpenSees.Elements.ElasticBeam2dWrapper(1, a, e, i, 1, 2, transftagCrdT, 0, 0, 0, 0);
            var ele21 = new OpenSees.Elements.ElasticBeam2dWrapper(21, a, e, i, 3, 6, transftagCrdT, 0, 0, 0, 0);
            var ele22 = new OpenSees.Elements.ElasticBeam2dWrapper(22, a, e, i, 6, 4, transftagCrdT, 0, 0, 0, 0);
            var ele3 = new OpenSees.Elements.ElasticBeam2dWrapper(3, a, e, i, 4, 5, transftagCrdT, 0, 0, 0, 0);
            theDomain.AddElement(ele1);
            theDomain.AddElement(ele21);
            theDomain.AddElement(ele22);
            theDomain.AddElement(ele3);

            //recorder Node -file data/node6disp.out -time -node 6 -dof 1 2 3 disp
            var savepath = System.Environment.CurrentDirectory + @"\opsnet_results_exp2\";
            if (!System.IO.Directory.Exists(savepath))
                System.IO.Directory.CreateDirectory(savepath);
            var dirIds2 = new OpenSees.IDWrapper(new int[] { 0, 1, 2 });
            var nodeTag = 6;
            var opsstream = new OpenSees.Handlers.DataFileStreamWrapper(savepath + $@"\OPSNET_Disp_{nodeTag}.txt");
            var recorder = new OpenSees.Recorders.NodeRecorderWrapper(dirIds2, new OpenSees.IDWrapper(new int[] { nodeTag }), 0, "disp", theDomain, opsstream);
            theDomain.AddRecorder(recorder);

            Console.WriteLine("model build!");


            // create linear pattern
            var theSeries = new OpenSees.Components.Timeseries.LinearSeriesWrapper();
            var theLoadPattern = new OpenSees.Components.LoadPatterns.LoadPatternWrapper(1);
            theLoadPattern.SetTimeSeries(theSeries);
            theDomain.AddLoadPattern(theLoadPattern);



            //$Wx mag of uniformily distributed ref load acting in direction along member length 
            //$Wy mag of uniformily distributed ref load acting in local y direction of element
            //$Wz mag of uniformily distributed ref load acting in local z direction of element
            //$Py mag of ref point load acting in direction along member length
            //$Py mag of ref point load acting in local y direction of element
            //$Pz mag of ref point load acting in local z direction of element
            //$xL location of point load relative to node I, prescribed as fraction of element length

            var wt = -5.0e4; // $Wy
            var wa = 0; // $Wx

            theDomain.AddElementLoad(new OpenSees.Components.Loads.Beam2dUniformLoadWrapper(1, wt, wa, 21), 1);
            theDomain.AddElementLoad(new OpenSees.Components.Loads.Beam2dUniformLoadWrapper(2, wt, wa, 22), 1);

            // for partial load use Beam2dPartialUniformLoadWrapper

            var theModel = new OpenSees.AnalysisModelWrapper();

            /*
                #define CURRENT_TANGENT 0
                #define INITIAL_TANGENT 1
                #define CURRENT_SECANT  2
                #define INITIAL_THEN_CURRENT_TANGENT  3
                #define NO_TANGENT  4
                #define SECOND_TANGENT 5
                #define HALL_TANGENT 6 
            */
            var theSolnAlgo = new OpenSees.Algorithms.ModifiedNewtonWrapper(0 /*CURRENT_TANGENT*/);

            var theIntegrator = new OpenSees.Integrators.Static.LoadControlWrapper(1, 1, 1, 1);
            var theHandler = new OpenSees.Handlers.TransformationConstraintHandlerWrapper();
            var theRCM = new OpenSees.GraphNumberers.RCMWrapper(false);
            var theNumberer = new OpenSees.Numberers.DOF_NumbererWrapper(theRCM);
            var theSolver = new OpenSees.Systems.Linears.BandGenLinLapackSolverWrapper();
            var theSOE = new OpenSees.Systems.Linears.BandGenLinSOEWrapper(theSolver);
            var test = new OpenSees.ConvergenceTests.CTestNormDispIncrWrapper(1e-8, 6, 2, 2, 1.0e10);
            var theAnalysis = new OpenSees.Analysis.StaticAnalysisWrapper(
                theDomain,
                theHandler,
                theNumberer,
                theModel,
                theSolnAlgo,
                theSOE,
                theIntegrator,
                test);

            theAnalysis.Analyze(1);


            // equal to loadConst -time 0.0
            theDomain.SetLoadConst(); // loadConst
            theDomain.SetCurrentTime(0); // -time 0.0
            theDomain.SetCommittedTime(0); // -time 0.0

            Console.ReadKey();
        }