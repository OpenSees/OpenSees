public static void Exp4()
        {
            Console.WriteLine("Exp4");
            var theDomain = new OpenSees.Components.DomainWrapper();

            var node1 = new OpenSees.Components.NodeWrapper(1, 3, 0, 0);
            var node2 = new OpenSees.Components.NodeWrapper(2, 3, 3.0, 0);
            theDomain.AddNode(new OpenSees.Components.NodeWrapper[] { node1, node2 });

            var sp1 = new OpenSees.Components.Constraints.SP_ConstraintWrapper(1, 0, 0.0, true);
            var sp2 = new OpenSees.Components.Constraints.SP_ConstraintWrapper(1, 1, 0.0, true);
            var sp3 = new OpenSees.Components.Constraints.SP_ConstraintWrapper(1, 2, 0.0, true);
            theDomain.AddSP_Constraint(new OpenSees.Components.Constraints.SP_ConstraintWrapper[] { sp1, sp2, sp3 });

            var E = 2.0e11;
            var steel01 = new OpenSees.Materials.Uniaxials.Steel01Wrapper(1, 2.354e8, E, 0.02);


            var repres = new OpenSees.Materials.Sections.Repres.FiberSectionReprWrapper(1,
                new OpenSees.Materials.Sections.Repres.PatchWrapper[]
                {
                        new OpenSees.Materials.Sections.Repres.QuadPatchWrapper(1,2,8,new OpenSees.MatrixWrapper(new double[,]{{ -0.15,0.125 },{-0.15,-0.125 },{ 0.15,-0.125 },{ 0.15,0.125, } })),
                },
                new OpenSees.Materials.Sections.Repres.ReinfLayerWrapper[]
                {
                });

            var section = new OpenSees.Materials.Sections.FiberSection2dWrapper(1, repres,
                new Dictionary<int, OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper>() {
                        {1,steel01 }
                });

            var goeTransCol = new OpenSees.Elements.CrdTransfs.LinearCrdTransf2dWrapper();
            var np = 5;
            var canti1 = new OpenSees.Elements.ForceBeamColumn2dWrapper(1, 1, 2, np, section, new OpenSees.Elements.BeamIntegrations.BeamIntegrationWrapper(OpenSees.Elements.BeamIntegrations.BeamIntegrationType.Lobatto), goeTransCol);
            theDomain.AddElement(new OpenSees.Elements.ElementWrapper[] { canti1 });

            {
                var opsstream = new OpenSees.Handlers.DataFileStreamWrapper(Environment.CurrentDirectory + $@"\exp4_disp.txt");
                var recorder = new OpenSees.Recorders.NodeRecorderWrapper(new OpenSees.IDWrapper(new int[] { 0, 1, 2 }), new OpenSees.IDWrapper(new int[] { 2 }), 0, "disp", theDomain, opsstream);
                theDomain.AddRecorder(recorder);
            }

            {
                var opsstream = new OpenSees.Handlers.DataFileStreamWrapper(Environment.CurrentDirectory + $@"\exp4_reactions.txt");
                var recorder = new OpenSees.Recorders.NodeRecorderWrapper(new OpenSees.IDWrapper(new int[] { 0, 1, 2 }), new OpenSees.IDWrapper(new int[] { 1 }), 0, "reaction", theDomain, opsstream);
                theDomain.AddRecorder(recorder);
            }


            var theSeries = new OpenSees.Components.Timeseries.LinearSeriesWrapper();
            var lpat = new OpenSees.Components.LoadPatterns.LoadPatternWrapper(3);
            lpat.SetTimeSeries(theSeries);
            theDomain.AddLoadPattern(lpat);

            var gravityLoadValues = new OpenSees.VectorWrapper(new double[] { 0, -100e3, 0 });
            var load = new OpenSees.Components.Loads.NodalLoadWrapper(1, 2, gravityLoadValues, false);
            theDomain.AddNodalLoad(new OpenSees.Components.Loads.NodalLoadWrapper[] { load }, 3);

            var theModel = new OpenSees.AnalysisModelWrapper();
            var theSolnAlgo = new OpenSees.Algorithms.NewtonRaphsonWrapper();
            var theIntegrator = new OpenSees.Integrators.Static.DisplacementControlWrapper(2, 1, -0.01, theDomain, 1, -0.01, -0.01);
            var theHandler = new OpenSees.Handlers.PlainHandlerWrapper();
            var theRCM = new OpenSees.GraphNumberers.RCMWrapper(false);
            var theNumberer = new OpenSees.Numberers.DOF_NumbererWrapper(theRCM);
            var theSolver = new OpenSees.Systems.Linears.BandSPDLinLapackSolverWrapper();
            var theSOE = new OpenSees.Systems.Linears.BandSPDLinSOEWrapper(theSolver);
            var theTest = new OpenSees.ConvergenceTests.CTestNormDispIncrWrapper(1e-8, 6, 2, 2, 1.0e10);
            var theAnalysis = new OpenSees.Analysis.StaticAnalysisWrapper(
                theDomain,
                theHandler,
                theNumberer,
                theModel,
                theSolnAlgo,
                theSOE,
                theIntegrator,
                theTest);

            theAnalysis.Analyze(50);


            node2.PrintSelf(0);
            
            Console.WriteLine("Exp4 completed");
        }