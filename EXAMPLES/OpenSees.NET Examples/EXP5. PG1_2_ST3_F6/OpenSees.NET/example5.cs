public static void Exp5()
        {
            var theDomain = new OpenSees.Components.DomainWrapper();
            #region parameters
            double length = 40.0;//set	length 									40.
            double width_inside_long = 40.0;//set width_inside_long						40.
            double width_inside_trans_Rec = 40.0;//set width_inside_trans_Rec					40.
            double width_Perimeter_long_right = 40.0;//set width_Perimeter_long_right				40.
            double width_Perimeter_long_left = 40.0;//set width_Perimeter_long_left				40.
            double width_Perimeter_trans_up = 40.0;//set width_Perimeter_trans_up				40.
            double width_Perimeter_trans_down = 40.0;//set width_Perimeter_trans_down				40.
            double tickness = 10.0;//set tickness								10.
            double Tension = 260000.0;//set Tension									[expr 260000.]
            double IP = 260000.0;//set IP										[expr 260000.]
            double OOP = 260000.0;//set OOP										[expr 260000.]
            double about_Tension = 260000.0;//set about_Tension							[expr 260000.]
            double about_IP = 260000.0;//set about_IP								[expr 260000.]
            double about_OOP = 260000.0;//set about_OOP								[expr 260000.]
            double E_Timber_Long = 13200.0; ;//set	E_Timber_Long							13200.0;
            double E_Timber_Long_mod = 15200.0; ;//set	E_Timber_Long_mod						15200.0;
            double E_Timber_Trans = 2200.0; ;//set	E_Timber_Trans							2200.0;
            double poisson_ratio = 0.04; ;//set	poisson_ratio							0.04;
            double G_Timber_Long = 240.0; ;//set	G_Timber_Long							240.0;
            double G_Timber_Trans = 220.0; ;//set	G_Timber_Trans							220.0;
            double Area_Perimeter_long_right = width_Perimeter_long_right * tickness;//set Area_Perimeter_long_right				[expr $width_Perimeter_long_right*$tickness]
            double Iz_Perimeter_long_right = (tickness * width_Perimeter_long_right * width_Perimeter_long_right * width_Perimeter_long_right) / (12.0);//set Iz_Perimeter_long_right					[expr ($tickness*$width_Perimeter_long_right*$width_Perimeter_long_right*$width_Perimeter_long_right)/(12.0)];
            double Iy_Perimeter_long_right = 1.0 * (width_Perimeter_long_right * tickness * tickness * tickness) / (12.0);//set Iy_Perimeter_long_right					[expr 1.0*($width_Perimeter_long_right*$tickness*$tickness*$tickness)/(12.0)];
            double Av_Perimeter_long_right = (5 * width_Perimeter_long_right * tickness) / (6);//set Av_Perimeter_long_right					[expr (5*$width_Perimeter_long_right*$tickness)/(6)];
            double Torsional_J_Perimeter_long_right = (width_Perimeter_long_right * tickness * (width_Perimeter_long_right * width_Perimeter_long_right + tickness * tickness)) / 12.0;//set Torsional_J_Perimeter_long_right		[expr ($width_Perimeter_long_right*$tickness*($width_Perimeter_long_right*$width_Perimeter_long_right+$tickness*$tickness))/12.0]
            double Iy_Perimeter_long_right_mod_down = 20.0 * (width_Perimeter_long_right * tickness * tickness * tickness) / (12.0);//set Iy_Perimeter_long_right_mod_down		[expr 20.0*($width_Perimeter_long_right*$tickness*$tickness*$tickness)/(12.0)];
            double Iy_Perimeter_long_right_mod_mid = 3.0 * (width_Perimeter_long_right * tickness * tickness * tickness) / (12.0);//set Iy_Perimeter_long_right_mod_mid			[expr 3.0*($width_Perimeter_long_right*$tickness*$tickness*$tickness)/(12.0)];
            double Iy_Perimeter_long_right_mod_up = 0.1 * (width_Perimeter_long_right * tickness * tickness * tickness) / (12.0);//set Iy_Perimeter_long_right_mod_up			[expr 0.1*($width_Perimeter_long_right*$tickness*$tickness*$tickness)/(12.0)];
            double Iz_Perimeter_long_right_mod_down = 48 * 100 * (tickness * width_Perimeter_long_right * width_Perimeter_long_right * width_Perimeter_long_right) / (12.0);//set Iz_Perimeter_long_right_mod_down		[expr 48*100*($tickness*$width_Perimeter_long_right*$width_Perimeter_long_right*$width_Perimeter_long_right)/(12.0)];
            double Iz_Perimeter_long_right_mod_mid = 48 * 100 * (tickness * width_Perimeter_long_right * width_Perimeter_long_right * width_Perimeter_long_right) / (12.0);//set Iz_Perimeter_long_right_mod_mid			[expr 48*100*($tickness*$width_Perimeter_long_right*$width_Perimeter_long_right*$width_Perimeter_long_right)/(12.0)];
            double Iz_Perimeter_long_right_mod_up = 48 * 100 * (tickness * width_Perimeter_long_right * width_Perimeter_long_right * width_Perimeter_long_right) / (12.0);//set Iz_Perimeter_long_right_mod_up			[expr 48*100*($tickness*$width_Perimeter_long_right*$width_Perimeter_long_right*$width_Perimeter_long_right)/(12.0)];
            double Area_Perimeter_long_left = width_Perimeter_long_left * tickness;//set Area_Perimeter_long_left				[expr $width_Perimeter_long_left*$tickness]
            double Iz_Perimeter_long_left = (tickness * width_Perimeter_long_left * width_Perimeter_long_left * width_Perimeter_long_left) / (12.0);//set Iz_Perimeter_long_left					[expr ($tickness*$width_Perimeter_long_left*$width_Perimeter_long_left*$width_Perimeter_long_left)/(12.0)];
            double Iy_Perimeter_long_left = (width_Perimeter_long_left * tickness * tickness * tickness) / (12.0);//set Iy_Perimeter_long_left					[expr ($width_Perimeter_long_left*$tickness*$tickness*$tickness)/(12.0)];
            double Av_Perimeter_long_left = (5 * width_Perimeter_long_left * tickness) / (6);//set Av_Perimeter_long_left					[expr (5*$width_Perimeter_long_left*$tickness)/(6)];
            double Torsional_J_Perimeter_long_left = (width_Perimeter_long_left * tickness * (width_Perimeter_long_left * width_Perimeter_long_left + tickness * tickness)) / 12.0;//set Torsional_J_Perimeter_long_left			[expr ($width_Perimeter_long_left*$tickness*($width_Perimeter_long_left*$width_Perimeter_long_left+$tickness*$tickness))/12.0]
            double Iz_Perimeter_long_left_mod_down = 0.25 * 0.01 * (tickness * width_Perimeter_long_left * width_Perimeter_long_left * width_Perimeter_long_left) / (12.0);//set Iz_Perimeter_long_left_mod_down			[expr 0.25*0.01*($tickness*$width_Perimeter_long_left*$width_Perimeter_long_left*$width_Perimeter_long_left)/(12.0)];
            double Iz_Perimeter_long_left_mod_mid = 0.25 * 0.1 * (tickness * width_Perimeter_long_left * width_Perimeter_long_left * width_Perimeter_long_left) / (12.0);//set Iz_Perimeter_long_left_mod_mid			[expr 0.25*0.1*($tickness*$width_Perimeter_long_left*$width_Perimeter_long_left*$width_Perimeter_long_left)/(12.0)];
            double Iz_Perimeter_long_left_mod_up = 0.25 * 0.1 * (tickness * width_Perimeter_long_left * width_Perimeter_long_left * width_Perimeter_long_left) / (12.0);//set Iz_Perimeter_long_left_mod_up			[expr 0.25*0.1*($tickness*$width_Perimeter_long_left*$width_Perimeter_long_left*$width_Perimeter_long_left)/(12.0)];
            double Area_Perimeter_trans_up = width_Perimeter_trans_up * tickness;//set Area_Perimeter_trans_up					[expr $width_Perimeter_trans_up*$tickness]
            double Iz_Perimeter_trans_up = (tickness * width_Perimeter_trans_up * width_Perimeter_trans_up * width_Perimeter_trans_up) / (12.0);//set Iz_Perimeter_trans_up					[expr ($tickness*$width_Perimeter_trans_up*$width_Perimeter_trans_up*$width_Perimeter_trans_up)/(12.0)];
            double Iy_Perimeter_trans_up = (width_Perimeter_trans_up * tickness * tickness * tickness) / (12.0);//set Iy_Perimeter_trans_up					[expr ($width_Perimeter_trans_up*$tickness*$tickness*$tickness)/(12.0)];
            double Av_Perimeter_trans_up = (5 * width_Perimeter_trans_up * tickness) / (6);//set Av_Perimeter_trans_up					[expr (5*$width_Perimeter_trans_up*$tickness)/(6)];
            double Torsional_J_Perimeter_trans_up = (width_Perimeter_trans_up * tickness * (width_Perimeter_trans_up * width_Perimeter_trans_up + tickness * tickness)) / 12.0;//set Torsional_J_Perimeter_trans_up			[expr ($width_Perimeter_trans_up*$tickness*($width_Perimeter_trans_up*$width_Perimeter_trans_up+$tickness*$tickness))/12.0]
            double Iy_Perimeter_trans_up_mod = 10.0 * (width_Perimeter_trans_up * tickness * tickness * tickness) / (12.0);//set Iy_Perimeter_trans_up_mod				[expr 10.0*($width_Perimeter_trans_up*$tickness*$tickness*$tickness)/(12.0)];
            double Iy_Perimeter_trans_up_mod_right = 0.5 * (width_Perimeter_trans_up * tickness * tickness * tickness) / (12.0);//set Iy_Perimeter_trans_up_mod_right			[expr 0.5*($width_Perimeter_trans_up*$tickness*$tickness*$tickness)/(12.0)];
            double Iz_Perimeter_trans_up_mod = 2.0 * (tickness * width_Perimeter_trans_up * width_Perimeter_trans_up * width_Perimeter_trans_up) / (12.0);//set Iz_Perimeter_trans_up_mod				[expr 2.*($tickness*$width_Perimeter_trans_up*$width_Perimeter_trans_up*$width_Perimeter_trans_up)/(12.0)];
            double Area_Perimeter_trans_down = width_Perimeter_trans_down * tickness;//set Area_Perimeter_trans_down				[expr $width_Perimeter_trans_down*$tickness]
            double Iz_Perimeter_trans_down = (tickness * width_Perimeter_trans_down * width_Perimeter_trans_down * width_Perimeter_trans_down) / (12.0);//set Iz_Perimeter_trans_down					[expr ($tickness*$width_Perimeter_trans_down*$width_Perimeter_trans_down*$width_Perimeter_trans_down)/(12.0)];
            double Iy_Perimeter_trans_down = (width_Perimeter_trans_down * tickness * tickness * tickness) / (12.0);//set Iy_Perimeter_trans_down					[expr ($width_Perimeter_trans_down*$tickness*$tickness*$tickness)/(12.0)];
            double Av_Perimeter_trans_down = (5 * width_Perimeter_trans_down * tickness) / (6);//set Av_Perimeter_trans_down					[expr (5*$width_Perimeter_trans_down*$tickness)/(6)];
            double Torsional_J_Perimeter_trans_down = (width_Perimeter_trans_down * tickness * (width_Perimeter_trans_down * width_Perimeter_trans_down + tickness * tickness)) / 12.0;//set Torsional_J_Perimeter_trans_down		[expr ($width_Perimeter_trans_down*$tickness*($width_Perimeter_trans_down*$width_Perimeter_trans_down+$tickness*$tickness))/12.0]
            double Iz_Perimeter_trans_down_mod = 0.20 * 0.775 * (tickness * width_Perimeter_trans_down * width_Perimeter_trans_down * width_Perimeter_trans_down) / (12.0);//set Iz_Perimeter_trans_down_mod				[expr 0.20*0.775*($tickness*$width_Perimeter_trans_down*$width_Perimeter_trans_down*$width_Perimeter_trans_down)/(12.0)];
            double Area_inside_long = width_inside_long * tickness;//set Area_inside_long						[expr $width_inside_long*$tickness]
            double Iz_inside_long = (tickness * width_inside_long * width_inside_long * width_inside_long) / (12.0);//set Iz_inside_long							[expr ($tickness*$width_inside_long*$width_inside_long*$width_inside_long)/(12.0)];
            double Iy_inside_long = 2.0 * (width_inside_long * tickness * tickness * tickness) / (12.0);//set Iy_inside_long							[expr 2.0*($width_inside_long*$tickness*$tickness*$tickness)/(12.0)];
            double Av_inside_long = (5 * width_inside_long * tickness) / (6);//set Av_inside_long							[expr (5*$width_inside_long*$tickness)/(6)];
            double Torsional_J_inside_long = (width_inside_long * tickness * (width_inside_long * width_inside_long + tickness * tickness)) / 12.0;//set Torsional_J_inside_long					[expr ($width_inside_long*$tickness*($width_inside_long*$width_inside_long+$tickness*$tickness))/12.0]
            double Iz_inside_long_mod = 1.0 * (tickness * width_inside_long * width_inside_long * width_inside_long) / (12.0);//set Iz_inside_long_mod						[expr 1.*($tickness*$width_inside_long*$width_inside_long*$width_inside_long)/(12.0)];
            double Area_inside_trans_Rec = width_inside_trans_Rec * tickness;//set Area_inside_trans_Rec					[expr $width_inside_trans_Rec*$tickness]
            double Iz_inside_trans_Rec = (tickness * width_inside_trans_Rec * width_inside_trans_Rec * width_inside_trans_Rec) / (12.0);//set Iz_inside_trans_Rec						[expr ($tickness*$width_inside_trans_Rec*$width_inside_trans_Rec*$width_inside_trans_Rec)/(12.0)];
            double Iy_inside_trans_Rec = 2.0 * (width_inside_trans_Rec * tickness * tickness * tickness) / (12.0);//set Iy_inside_trans_Rec						[expr 2.0*($width_inside_trans_Rec*$tickness*$tickness*$tickness)/(12.0)];
            double Av_inside_trans_Rec = (5 * width_inside_trans_Rec * tickness) / (6);//set Av_inside_trans_Rec						[expr (5*$width_inside_trans_Rec*$tickness)/(6)];
            double Torsional_J_inside_trans_Rec = (width_inside_trans_Rec * tickness * (width_inside_trans_Rec * width_inside_trans_Rec + tickness * tickness)) / 12.0;//set Torsional_J_inside_trans_Rec			[expr ($width_inside_trans_Rec*$tickness*($width_inside_trans_Rec*$width_inside_trans_Rec+$tickness*$tickness))/12.0]
            double Iz_inside_trans_Rec_mod = 1.0 * (tickness * width_inside_trans_Rec * width_inside_trans_Rec * width_inside_trans_Rec) / (12.0);//set Iz_inside_trans_Rec_mod					[expr 1.*($tickness*$width_inside_trans_Rec*$width_inside_trans_Rec*$width_inside_trans_Rec)/(12.0)];

            #endregion

            #region materials
            int rigidmattag = 1;
            int semi_rigidmattag = 101;
            int Freemattag = 2;
            int timbermattag_long = 30;
            int timbermattag_trans = 31;
            int G_Timber_Tag_Long = 32;
            int G_Timber_Tag_Trans = 33;
            int timbermattag_Ortho_Long = 21;
            int timbermattag_Ortho_Trans = 22;
            int Stiffness_Tension = 10;
            int Stiffness_IP = 11;
            int Stiffness_OOP = 12;
            int Stiffness_about_Tension = 13;
            int Stiffness_about_IP = 14;
            int Stiffness_about_OOP = 15;
            int Trans_Shear_Timoshenko_Eq = 70;
            int Long_Shear_Timoshenko_Eq = 71;

            var rigidmattagMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(rigidmattag, 1.0e+15, 0);//uniaxialMaterial	Elastic	$rigidmattag  						[expr 1.0e+15]
            var semi_rigidmattagMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(semi_rigidmattag, 1.0e+4, 0);//uniaxialMaterial	Elastic	$semi_rigidmattag					[expr 1.0e+4]
            var FreemattagMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(Freemattag, 1.0e-1, 0);//uniaxialMaterial	Elastic	$Freemattag  						[expr 1.0e-1]
            var timbermattag_longMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(timbermattag_long, E_Timber_Long, 0);//uniaxialMaterial	Elastic	$timbermattag_long					$E_Timber_Long
            var timbermattag_transMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(timbermattag_trans, E_Timber_Trans, 0);//uniaxialMaterial	Elastic	$timbermattag_trans					$E_Timber_Trans
            var G_Timber_Tag_LongMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(G_Timber_Tag_Long, G_Timber_Long, 0);//uniaxialMaterial	Elastic	$G_Timber_Tag_Long					$G_Timber_Long
            var G_Timber_Tag_TransMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(G_Timber_Tag_Trans, G_Timber_Trans, 0);//uniaxialMaterial	Elastic	$G_Timber_Tag_Trans					$G_Timber_Trans
            var Stiffness_TensionMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(Stiffness_Tension, Tension, 0);//uniaxialMaterial	Elastic	$Stiffness_Tension					[expr $Tension]
            var Stiffness_IPMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(Stiffness_IP, IP, 0);//uniaxialMaterial	Elastic	$Stiffness_IP						[expr $IP]
            var Stiffness_OOPMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(Stiffness_OOP, OOP, 0);//uniaxialMaterial	Elastic	$Stiffness_OOP						[expr $OOP]
            var Stiffness_about_TensionMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(Stiffness_about_Tension, about_Tension, 0);//uniaxialMaterial	Elastic	$Stiffness_about_Tension			[expr $about_Tension]
            var Stiffness_about_IPMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(Stiffness_about_IP, about_IP, 0);//uniaxialMaterial	Elastic	$Stiffness_about_IP					[expr $about_IP]
            var Stiffness_about_OOPMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(Stiffness_about_OOP, about_OOP, 0);//uniaxialMaterial	Elastic	$Stiffness_about_OOP				[expr $about_OOP]
            var Trans_Shear_Timoshenko_EqMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(Trans_Shear_Timoshenko_Eq, 34787.878, 0);//uniaxialMaterial	Elastic	$Trans_Shear_Timoshenko_Eq	[expr 34787.878]
            var Long_Shear_Timoshenko_EqMat = new OpenSees.Materials.Uniaxials.ElasticMaterialWrapper(Long_Shear_Timoshenko_Eq, 16081.971, 0);//uniaxialMaterial	Elastic	$Long_Shear_Timoshenko_Eq	[expr 16081.971]

            #endregion
            #region nodes
            var node401 = new OpenSees.Components.NodeWrapper(401, 6, 880.0, 0.0, 1400.0);//node	401		880.		0.0		1400.
            theDomain.AddNode(node401);

            var node402 = new OpenSees.Components.NodeWrapper(402, 6, 697.394, 0.0, 1400.0);//node	402		697.394		0.0		1400.
            theDomain.AddNode(node402);

            var node403 = new OpenSees.Components.NodeWrapper(403, 6, 684.894, 0.0, 1400.0);//node	403		684.894		0.0		1400.
            theDomain.AddNode(node403);

            var node404 = new OpenSees.Components.NodeWrapper(404, 6, 672.394, 0.0, 1400.0);//node	404		672.394		0.0		1400.
            theDomain.AddNode(node404);

            var node405 = new OpenSees.Components.NodeWrapper(405, 6, 660.0, 0.0, 1400.0);//node	405		660.		0.0		1400.
            theDomain.AddNode(node405);

            var node406 = new OpenSees.Components.NodeWrapper(406, 6, 659.894, 0.0, 1400.0);//node	406		659.894		0.0		1400.
            theDomain.AddNode(node406);

            var node407 = new OpenSees.Components.NodeWrapper(407, 6, 647.5, 0.0, 1400.0);//node	407		647.5		0.0		1400.
            theDomain.AddNode(node407);

            var node408 = new OpenSees.Components.NodeWrapper(408, 6, 440.0, 0.0, 1400.0);//node	408		440.		0.0		1400.
            theDomain.AddNode(node408);

            var node409 = new OpenSees.Components.NodeWrapper(409, 6, 280.207, 0.0, 1400.0);//node	409		280.207		0.0		1400.
            theDomain.AddNode(node409);

            var node410 = new OpenSees.Components.NodeWrapper(410, 6, 268.957, 0.0, 1400.0);//node	410		268.957		0.0		1400.
            theDomain.AddNode(node410);

            var node411 = new OpenSees.Components.NodeWrapper(411, 6, 267.707, 0.0, 1400.0);//node	411		267.707		0.0		1400.
            theDomain.AddNode(node411);

            var node412 = new OpenSees.Components.NodeWrapper(412, 6, 255.207, 0.0, 1400.0);//node	412		255.207		0.0		1400.
            theDomain.AddNode(node412);

            var node413 = new OpenSees.Components.NodeWrapper(413, 6, 242.707, 0.0, 1400.0);//node	413		242.707		0.0		1400.
            theDomain.AddNode(node413);

            var node414 = new OpenSees.Components.NodeWrapper(414, 6, 230.207, 0.0, 1400.0);//node	414		230.207		0.0		1400.
            theDomain.AddNode(node414);

            var node415 = new OpenSees.Components.NodeWrapper(415, 6, 220.0, 0.0, 1400.0);//node	415		220.		0.0		1400.
            theDomain.AddNode(node415);

            var node4020 = new OpenSees.Components.NodeWrapper(4020, 6, 697.394, 0.0, 1440.0);//node	4020	697.394		0.0		1440.
            theDomain.AddNode(node4020);

            var node4030 = new OpenSees.Components.NodeWrapper(4030, 6, 684.894, 0.0, 1440.0);//node	4030	684.894		0.0		1440.
            theDomain.AddNode(node4030);

            var node4040 = new OpenSees.Components.NodeWrapper(4040, 6, 672.394, 0.0, 1440.0);//node	4040	672.394		0.0		1440.
            theDomain.AddNode(node4040);

            var node4060 = new OpenSees.Components.NodeWrapper(4060, 6, 659.894, 0.0, 1440.0);//node	4060	659.894		0.0		1440.
            theDomain.AddNode(node4060);

            var node4070 = new OpenSees.Components.NodeWrapper(4070, 6, 647.5, 0.0, 1440.0);//node	4070	647.5		0.0		1440.
            theDomain.AddNode(node4070);

            var node4090 = new OpenSees.Components.NodeWrapper(4090, 6, 280.207, 0.0, 1440.0);//node	4090	280.207		0.0		1440.
            theDomain.AddNode(node4090);

            var node4110 = new OpenSees.Components.NodeWrapper(4110, 6, 267.707, 0.0, 1440.0);//node	4110	267.707		0.0		1440.
            theDomain.AddNode(node4110);

            var node4130 = new OpenSees.Components.NodeWrapper(4130, 6, 242.707, 0.0, 1440.0);//node	4130	242.707		0.0		1440.
            theDomain.AddNode(node4130);

            var node4140 = new OpenSees.Components.NodeWrapper(4140, 6, 230.207, 0.0, 1440.0);//node	4140	230.207		0.0		1440.
            theDomain.AddNode(node4140);

            var node4120 = new OpenSees.Components.NodeWrapper(4120, 6, 255.207, 0.0, 1440.0);//node	4120	255.207		0.0		1440.
            theDomain.AddNode(node4120);

            var node301 = new OpenSees.Components.NodeWrapper(301, 6, 880.0, 0.0, 0.0);//node	301		880.		0.0		0.0
            theDomain.AddNode(node301);

            var node302 = new OpenSees.Components.NodeWrapper(302, 6, 863.75, 0.0, 0.0);//node	302		863.75		0.0		0.
            theDomain.AddNode(node302);

            var node303 = new OpenSees.Components.NodeWrapper(303, 6, 852.5, 0.0, 0.0);//node	303		852.5		0.0		0.
            theDomain.AddNode(node303);

            var node304 = new OpenSees.Components.NodeWrapper(304, 6, 851.25, 0.0, 0.0);//node	304		851.25		0.0		0.
            theDomain.AddNode(node304);

            var node305 = new OpenSees.Components.NodeWrapper(305, 6, 838.75, 0.0, 0.0);//node	305		838.75		0.0		0.
            theDomain.AddNode(node305);

            var node306 = new OpenSees.Components.NodeWrapper(306, 6, 826.25, 0.0, 0.0);//node	306		826.25		0.0		0.
            theDomain.AddNode(node306);

            var node307 = new OpenSees.Components.NodeWrapper(307, 6, 813.75, 0.0, 0.0);//node	307		813.75		0.0		0.
            theDomain.AddNode(node307);

            var node308 = new OpenSees.Components.NodeWrapper(308, 6, 660.0, 0.0, 0.0);//node	308		660.		0.0		0.0
            theDomain.AddNode(node308);

            var node309 = new OpenSees.Components.NodeWrapper(309, 6, 451.25, 0.0, 0.0);//node	309		451.25		0.0		0.
            theDomain.AddNode(node309);

            var node310 = new OpenSees.Components.NodeWrapper(310, 6, 440.0, 0.0, 0.0);//node	310		440.		0.0		0.0
            theDomain.AddNode(node310);

            var node311 = new OpenSees.Components.NodeWrapper(311, 6, 438.75, 0.0, 0.0);//node	311		438.75		0.0		0.
            theDomain.AddNode(node311);

            var node312 = new OpenSees.Components.NodeWrapper(312, 6, 426.25, 0.0, 0.0);//node	312		426.25		0.0		0.
            theDomain.AddNode(node312);

            var node313 = new OpenSees.Components.NodeWrapper(313, 6, 413.75, 0.0, 0.0);//node	313		413.75		0.0		0.
            theDomain.AddNode(node313);

            var node314 = new OpenSees.Components.NodeWrapper(314, 6, 401.25, 0.0, 0.0);//node	314		401.25		0.0		0.
            theDomain.AddNode(node314);

            var node315 = new OpenSees.Components.NodeWrapper(315, 6, 220.0, 0.0, 0.0);//node	315		220.		0.0		0.0
            theDomain.AddNode(node315);

            var node3020 = new OpenSees.Components.NodeWrapper(3020, 6, 863.75, 0.0, -40.0);//node	3020	863.75		0.0		-40.
            theDomain.AddNode(node3020);

            var node3040 = new OpenSees.Components.NodeWrapper(3040, 6, 851.25, 0.0, -40.0);//node	3040	851.25		0.0		-40.
            theDomain.AddNode(node3040);

            var node3050 = new OpenSees.Components.NodeWrapper(3050, 6, 838.75, 0.0, -40.0);//node	3050	838.75		0.0		-40.
            theDomain.AddNode(node3050);

            var node3060 = new OpenSees.Components.NodeWrapper(3060, 6, 826.25, 0.0, -40.0);//node	3060	826.25		0.0		-40.
            theDomain.AddNode(node3060);

            var node3070 = new OpenSees.Components.NodeWrapper(3070, 6, 813.75, 0.0, -40.0);//node	3070	813.75		0.0		-40.
            theDomain.AddNode(node3070);

            var node3090 = new OpenSees.Components.NodeWrapper(3090, 6, 451.25, 0.0, -40.0);//node	3090	451.25		0.0		-40.
            theDomain.AddNode(node3090);

            var node3110 = new OpenSees.Components.NodeWrapper(3110, 6, 438.75, 0.0, -40.0);//node	3110	438.75		0.0		-40.
            theDomain.AddNode(node3110);

            var node3120 = new OpenSees.Components.NodeWrapper(3120, 6, 426.25, 0.0, -40.0);//node	3120	426.25		0.0		-40.
            theDomain.AddNode(node3120);

            var node3130 = new OpenSees.Components.NodeWrapper(3130, 6, 413.75, 0.0, -40.0);//node	3130	413.75		0.0		-40.
            theDomain.AddNode(node3130);

            var node3140 = new OpenSees.Components.NodeWrapper(3140, 6, 401.25, 0.0, -40.0);//node	3140	401.25		0.0		-40.
            theDomain.AddNode(node3140);

            var node200 = new OpenSees.Components.NodeWrapper(200, 6, 1100.0, 0.0, 0.0);//node	200		1100.		0.0		0.
            theDomain.AddNode(node200);

            var node201 = new OpenSees.Components.NodeWrapper(201, 6, 1100.0, 0.0, 280.0);//node	201		1100.		0.0		280.
            theDomain.AddNode(node201);

            var node202 = new OpenSees.Components.NodeWrapper(202, 6, 1100.0, 0.0, 330.0);//node	202		1100.		0.0		330.
            theDomain.AddNode(node202);

            var node203 = new OpenSees.Components.NodeWrapper(203, 6, 1100.0, 0.0, 342.5);//node	203		1100.		0.0		342.5
            theDomain.AddNode(node203);

            var node204 = new OpenSees.Components.NodeWrapper(204, 6, 1100.0, 0.0, 355.0);//node	204		1100.		0.0		355.
            theDomain.AddNode(node204);

            var node205 = new OpenSees.Components.NodeWrapper(205, 6, 1100.0, 0.0, 367.5);//node	205		1100.		0.0		367.5
            theDomain.AddNode(node205);

            var node206 = new OpenSees.Components.NodeWrapper(206, 6, 1100.0, 0.0, 380.0);//node	206		1100.		0.0		380.
            theDomain.AddNode(node206);

            var node207 = new OpenSees.Components.NodeWrapper(207, 6, 1100.0, 0.0, 560.0);//node	207		1100.		0.0		560.
            theDomain.AddNode(node207);

            var node208 = new OpenSees.Components.NodeWrapper(208, 6, 1100.0, 0.0, 740.0);//node	208		1100.		0.0		740.
            theDomain.AddNode(node208);

            var node209 = new OpenSees.Components.NodeWrapper(209, 6, 1100.0, 0.0, 752.5);//node	209		1100.		0.0		752.5
            theDomain.AddNode(node209);

            var node210 = new OpenSees.Components.NodeWrapper(210, 6, 1100.0, 0.0, 765.0);//node	210		1100.		0.0		765.
            theDomain.AddNode(node210);

            var node211 = new OpenSees.Components.NodeWrapper(211, 6, 1100.0, 0.0, 777.5);//node	211		1100.		0.0		777.5
            theDomain.AddNode(node211);

            var node212 = new OpenSees.Components.NodeWrapper(212, 6, 1100.0, 0.0, 790.0);//node	212		1100.		0.0		790.
            theDomain.AddNode(node212);

            var node213 = new OpenSees.Components.NodeWrapper(213, 6, 1100.0, 0.0, 840.0);//node	213		1100.		0.0		840.
            theDomain.AddNode(node213);

            var node214 = new OpenSees.Components.NodeWrapper(214, 6, 1100.0, 0.0, 1120.0);//node	214		1100.		0.0		1120.
            theDomain.AddNode(node214);

            var node215 = new OpenSees.Components.NodeWrapper(215, 6, 1100.0, 0.0, 1130.0);//node	215		1100.		0.0		1130.
            theDomain.AddNode(node215);

            var node216 = new OpenSees.Components.NodeWrapper(216, 6, 1100.0, 0.0, 1142.5);//node	216		1100.		0.0		1142.5
            theDomain.AddNode(node216);

            var node217 = new OpenSees.Components.NodeWrapper(217, 6, 1100.0, 0.0, 1155.0);//node	217		1100.		0.0		1155.
            theDomain.AddNode(node217);

            var node218 = new OpenSees.Components.NodeWrapper(218, 6, 1100.0, 0.0, 1167.5);//node	218		1100.		0.0		1167.5
            theDomain.AddNode(node218);

            var node219 = new OpenSees.Components.NodeWrapper(219, 6, 1100.0, 0.0, 1180.0);//node	219		1100.		0.0		1180.
            theDomain.AddNode(node219);

            var node220 = new OpenSees.Components.NodeWrapper(220, 6, 1100.0, 0.0, 1400.0);//node	220		1100.		0.0		1400.
            theDomain.AddNode(node220);

            var node2020 = new OpenSees.Components.NodeWrapper(2020, 6, 1140.0, 0.0, 330.0);//node	2020	1140.		0.0		330.
            theDomain.AddNode(node2020);

            var node2030 = new OpenSees.Components.NodeWrapper(2030, 6, 1140.0, 0.0, 342.5);//node	2030	1140.		0.0		342.5
            theDomain.AddNode(node2030);

            var node2040 = new OpenSees.Components.NodeWrapper(2040, 6, 1140, 0.0, 355.0);//node	2040	1140		0.0		355.
            theDomain.AddNode(node2040);

            var node2050 = new OpenSees.Components.NodeWrapper(2050, 6, 1140.0, 0.0, 367.5);//node	2050	1140.		0.0		367.5
            theDomain.AddNode(node2050);

            var node2060 = new OpenSees.Components.NodeWrapper(2060, 6, 1140.0, 0.0, 380.0);//node	2060	1140.		0.0		380.
            theDomain.AddNode(node2060);

            var node2080 = new OpenSees.Components.NodeWrapper(2080, 6, 1140.0, 0.0, 740.0);//node	2080	1140.		0.0		740.
            theDomain.AddNode(node2080);

            var node2090 = new OpenSees.Components.NodeWrapper(2090, 6, 1140.0, 0.0, 752.5);//node	2090	1140.		0.0		752.5
            theDomain.AddNode(node2090);

            var node2100 = new OpenSees.Components.NodeWrapper(2100, 6, 1140.0, 0.0, 765.0);//node	2100	1140.		0.0		765.
            theDomain.AddNode(node2100);

            var node2110 = new OpenSees.Components.NodeWrapper(2110, 6, 1140.0, 0.0, 777.5);//node	2110	1140.		0.0		777.5
            theDomain.AddNode(node2110);

            var node2120 = new OpenSees.Components.NodeWrapper(2120, 6, 1140.0, 0.0, 790.0);//node	2120	1140.		0.0		790.
            theDomain.AddNode(node2120);

            var node2150 = new OpenSees.Components.NodeWrapper(2150, 6, 1140.0, 0.0, 1130.0);//node	2150	1140.		0.0		1130.
            theDomain.AddNode(node2150);

            var node2160 = new OpenSees.Components.NodeWrapper(2160, 6, 1140.0, 0.0, 1142.5);//node	2160	1140.		0.0		1142.5
            theDomain.AddNode(node2160);

            var node2170 = new OpenSees.Components.NodeWrapper(2170, 6, 1140.0, 0.0, 1155.0);//node	2170	1140.		0.0		1155.
            theDomain.AddNode(node2170);

            var node2180 = new OpenSees.Components.NodeWrapper(2180, 6, 1140.0, 0.0, 1167.5);//node	2180	1140.		0.0		1167.5
            theDomain.AddNode(node2180);

            var node2190 = new OpenSees.Components.NodeWrapper(2190, 6, 1140.0, 0.0, 1180.0);//node	2190	1140.		0.0		1180.
            theDomain.AddNode(node2190);

            var node100 = new OpenSees.Components.NodeWrapper(100, 6, 0.0, 0.0, 0.0);//node	100		0.			0.0		0.
            theDomain.AddNode(node100);

            var node101 = new OpenSees.Components.NodeWrapper(101, 6, 0.0, 0.0, 220.0);//node	101		0.			0.0		220.
            theDomain.AddNode(node101);

            var node102 = new OpenSees.Components.NodeWrapper(102, 6, 0.0, 0.0, 232.5);//node	102		0.			0.0		232.5
            theDomain.AddNode(node102);

            var node103 = new OpenSees.Components.NodeWrapper(103, 6, 0.0, 0.0, 245.0);//node	103		0.			0.0		245.
            theDomain.AddNode(node103);

            var node104 = new OpenSees.Components.NodeWrapper(104, 6, 0.0, 0.0, 257.5);//node	104		0.			0.0		257.5
            theDomain.AddNode(node104);

            var node105 = new OpenSees.Components.NodeWrapper(105, 6, 0.0, 0.0, 270.0);//node	105		0.			0.0		270.
            theDomain.AddNode(node105);

            var node106 = new OpenSees.Components.NodeWrapper(106, 6, 0.0, 0.0, 280.0);//node	106		0.			0.0		280.
            theDomain.AddNode(node106);

            var node107 = new OpenSees.Components.NodeWrapper(107, 6, 0.0, 0.0, 560.0);//node	107		0.			0.0		560.
            theDomain.AddNode(node107);

            var node108 = new OpenSees.Components.NodeWrapper(108, 6, 0.0, 0.0, 630.0);//node	108		0.			0.0		630.
            theDomain.AddNode(node108);

            var node109 = new OpenSees.Components.NodeWrapper(109, 6, 0.0, 0.0, 642.5);//node	109		0.			0.0		642.5
            theDomain.AddNode(node109);

            var node110 = new OpenSees.Components.NodeWrapper(110, 6, 0.0, 0.0, 655.0);//node	110		0.			0.0		655.
            theDomain.AddNode(node110);

            var node111 = new OpenSees.Components.NodeWrapper(111, 6, 0.0, 0.0, 667.5);//node	111		0.			0.0		667.5
            theDomain.AddNode(node111);

            var node112 = new OpenSees.Components.NodeWrapper(112, 6, 0.0, 0.0, 680.0);//node	112		0.			0.0		680.
            theDomain.AddNode(node112);

            var node113 = new OpenSees.Components.NodeWrapper(113, 6, 0.0, 0.0, 840.0);//node	113		0.			0.0		840.
            theDomain.AddNode(node113);

            var node114 = new OpenSees.Components.NodeWrapper(114, 6, 0.0, 0.0, 1020.0);//node	114		0.			0.0		1020.
            theDomain.AddNode(node114);

            var node115 = new OpenSees.Components.NodeWrapper(115, 6, 0.0, 0.0, 1032.5);//node	115		0.			0.0		1032.5
            theDomain.AddNode(node115);

            var node116 = new OpenSees.Components.NodeWrapper(116, 6, 0.0, 0.0, 1045.0);//node	116		0.			0.0		1045.
            theDomain.AddNode(node116);

            var node117 = new OpenSees.Components.NodeWrapper(117, 6, 0.0, 0.0, 1057.5);//node	117		0.			0.0		1057.5
            theDomain.AddNode(node117);

            var node118 = new OpenSees.Components.NodeWrapper(118, 6, 0.0, 0.0, 1070.0);//node	118		0.			0.0		1070.
            theDomain.AddNode(node118);

            var node119 = new OpenSees.Components.NodeWrapper(119, 6, 0.0, 0.0, 1120.0);//node	119		0.			0.0		1120.
            theDomain.AddNode(node119);

            var node120 = new OpenSees.Components.NodeWrapper(120, 6, 0.0, 0.0, 1400.0);//node	120		0.			0.0		1400.
            theDomain.AddNode(node120);

            var node1010 = new OpenSees.Components.NodeWrapper(1010, 6, -40.0, 0.0, 220.0);//node	1010	-40.		0.0		220.
            theDomain.AddNode(node1010);

            var node1020 = new OpenSees.Components.NodeWrapper(1020, 6, -40.0, 0.0, 232.5);//node	1020	-40.		0.0		232.5
            theDomain.AddNode(node1020);

            var node1030 = new OpenSees.Components.NodeWrapper(1030, 6, -40.0, 0.0, 245.0);//node	1030	-40.		0.0		245.
            theDomain.AddNode(node1030);

            var node1040 = new OpenSees.Components.NodeWrapper(1040, 6, -40.0, 0.0, 257.5);//node	1040	-40.		0.0		257.5
            theDomain.AddNode(node1040);

            var node1050 = new OpenSees.Components.NodeWrapper(1050, 6, -40.0, 0.0, 270.0);//node	1050	-40.		0.0		270.
            theDomain.AddNode(node1050);

            var node1080 = new OpenSees.Components.NodeWrapper(1080, 6, -40.0, 0.0, 630.0);//node	1080	-40.		0.0		630.
            theDomain.AddNode(node1080);

            var node1090 = new OpenSees.Components.NodeWrapper(1090, 6, -40.0, 0.0, 642.5);//node	1090	-40.		0.0		642.5
            theDomain.AddNode(node1090);

            var node1100 = new OpenSees.Components.NodeWrapper(1100, 6, -40.0, 0.0, 655.0);//node	1100	-40.		0.0		655.
            theDomain.AddNode(node1100);

            var node1110 = new OpenSees.Components.NodeWrapper(1110, 6, -40.0, 0.0, 667.5);//node	1110	-40.		0.0		667.5
            theDomain.AddNode(node1110);

            var node1120 = new OpenSees.Components.NodeWrapper(1120, 6, -40.0, 0.0, 680.0);//node	1120	-40.		0.0		680.
            theDomain.AddNode(node1120);

            var node1140 = new OpenSees.Components.NodeWrapper(1140, 6, -40.0, 0.0, 1020.0);//node	1140	-40.		0.0		1020.
            theDomain.AddNode(node1140);

            var node1150 = new OpenSees.Components.NodeWrapper(1150, 6, -40.0, 0.0, 1032.5);//node	1150	-40.		0.0		1032.5
            theDomain.AddNode(node1150);

            var node1160 = new OpenSees.Components.NodeWrapper(1160, 6, -40.0, 0.0, 1045.0);//node	1160	-40.		0.0		1045.
            theDomain.AddNode(node1160);

            var node1170 = new OpenSees.Components.NodeWrapper(1170, 6, -40.0, 0.0, 1057.5);//node	1170	-40.		0.0		1057.5
            theDomain.AddNode(node1170);

            var node1180 = new OpenSees.Components.NodeWrapper(1180, 6, -40.0, 0.0, 1070.0);//node	1180	-40.		0.0		1070.
            theDomain.AddNode(node1180);

            var node10001 = new OpenSees.Components.NodeWrapper(10001, 6, 0.0, 0.0, 2000 + 0.0);//node	10001		0.			0.0		[expr 2000+0.]
            theDomain.AddNode(node10001);

            var node10101 = new OpenSees.Components.NodeWrapper(10101, 6, 0.0, 0.0, 2000 + 220.0);//node	10101		0.			0.0		[expr 2000+220.]
            theDomain.AddNode(node10101);

            var node10201 = new OpenSees.Components.NodeWrapper(10201, 6, 0.0, 0.0, 2000 + 232.5);//node	10201		0.			0.0		[expr 2000+232.5]
            theDomain.AddNode(node10201);

            var node10301 = new OpenSees.Components.NodeWrapper(10301, 6, 0.0, 0.0, 2000 + 245.0);//node	10301		0.			0.0		[expr 2000+245.]
            theDomain.AddNode(node10301);

            var node10401 = new OpenSees.Components.NodeWrapper(10401, 6, 0.0, 0.0, 2000 + 257.5);//node	10401		0.			0.0		[expr 2000+257.5]
            theDomain.AddNode(node10401);

            var node10501 = new OpenSees.Components.NodeWrapper(10501, 6, 0.0, 0.0, 2000 + 270.0);//node	10501		0.			0.0		[expr 2000+270.]
            theDomain.AddNode(node10501);

            var node10601 = new OpenSees.Components.NodeWrapper(10601, 6, 0.0, 0.0, 2000 + 280.0);//node	10601		0.			0.0		[expr 2000+280.]
            theDomain.AddNode(node10601);

            var node10701 = new OpenSees.Components.NodeWrapper(10701, 6, 0.0, 0.0, 2000 + 560.0);//node	10701		0.			0.0		[expr 2000+560.]
            theDomain.AddNode(node10701);

            var node10801 = new OpenSees.Components.NodeWrapper(10801, 6, 0.0, 0.0, 2000 + 630.0);//node	10801		0.			0.0		[expr 2000+630.]
            theDomain.AddNode(node10801);

            var node10901 = new OpenSees.Components.NodeWrapper(10901, 6, 0.0, 0.0, 2000 + 642.5);//node	10901		0.			0.0		[expr 2000+642.5]
            theDomain.AddNode(node10901);

            var node11001 = new OpenSees.Components.NodeWrapper(11001, 6, 0.0, 0.0, 2000 + 655.0);//node	11001		0.			0.0		[expr 2000+655.]
            theDomain.AddNode(node11001);

            var node11101 = new OpenSees.Components.NodeWrapper(11101, 6, 0.0, 0.0, 2000 + 667.5);//node	11101		0.			0.0		[expr 2000+667.5]
            theDomain.AddNode(node11101);

            var node11201 = new OpenSees.Components.NodeWrapper(11201, 6, 0.0, 0.0, 2000 + 680.0);//node	11201		0.			0.0		[expr 2000+680.]
            theDomain.AddNode(node11201);

            var node11301 = new OpenSees.Components.NodeWrapper(11301, 6, 0.0, 0.0, 2000 + 840.0);//node	11301		0.			0.0		[expr 2000+840.]
            theDomain.AddNode(node11301);

            var node11401 = new OpenSees.Components.NodeWrapper(11401, 6, 0.0, 0.0, 2000 + 1020.0);//node	11401		0.			0.0		[expr 2000+1020.]
            theDomain.AddNode(node11401);

            var node11501 = new OpenSees.Components.NodeWrapper(11501, 6, 0.0, 0.0, 2000 + 1032.5);//node	11501		0.			0.0		[expr 2000+1032.5]
            theDomain.AddNode(node11501);

            var node11601 = new OpenSees.Components.NodeWrapper(11601, 6, 0.0, 0.0, 2000 + 1045.0);//node	11601		0.			0.0		[expr 2000+1045.]
            theDomain.AddNode(node11601);

            var node11701 = new OpenSees.Components.NodeWrapper(11701, 6, 0.0, 0.0, 2000 + 1057.5);//node	11701		0.			0.0		[expr 2000+1057.5]
            theDomain.AddNode(node11701);

            var node11801 = new OpenSees.Components.NodeWrapper(11801, 6, 0.0, 0.0, 2000 + 1070.0);//node	11801		0.			0.0		[expr 2000+1070.]
            theDomain.AddNode(node11801);

            var node11901 = new OpenSees.Components.NodeWrapper(11901, 6, 0.0, 0.0, 2000 + 1120.0);//node	11901		0.			0.0		[expr 2000+1120.]
            theDomain.AddNode(node11901);

            var node12001 = new OpenSees.Components.NodeWrapper(12001, 6, 0.0, 0.0, 2000 + 1400.0);//node	12001		0.			0.0		[expr 2000+1400.]
            theDomain.AddNode(node12001);

            var node20001 = new OpenSees.Components.NodeWrapper(20001, 6, 1100.0, 0.0, 2000 + 0.0);//node	20001		1100.		0.0		[expr 2000+0.]
            theDomain.AddNode(node20001);

            var node20101 = new OpenSees.Components.NodeWrapper(20101, 6, 1100.0, 0.0, 2000 + 280.0);//node	20101		1100.		0.0		[expr 2000+280.]
            theDomain.AddNode(node20101);

            var node20201 = new OpenSees.Components.NodeWrapper(20201, 6, 1100.0, 0.0, 2000 + 330.0);//node	20201		1100.		0.0		[expr 2000+330.]
            theDomain.AddNode(node20201);

            var node20301 = new OpenSees.Components.NodeWrapper(20301, 6, 1100.0, 0.0, 2000 + 342.5);//node	20301		1100.		0.0		[expr 2000+342.5]
            theDomain.AddNode(node20301);

            var node20401 = new OpenSees.Components.NodeWrapper(20401, 6, 1100.0, 0.0, 2000 + 355.0);//node	20401		1100.		0.0		[expr 2000+355.]
            theDomain.AddNode(node20401);

            var node20501 = new OpenSees.Components.NodeWrapper(20501, 6, 1100.0, 0.0, 2000 + 367.5);//node	20501		1100.		0.0		[expr 2000+367.5]
            theDomain.AddNode(node20501);

            var node20601 = new OpenSees.Components.NodeWrapper(20601, 6, 1100.0, 0.0, 2000 + 380.0);//node	20601		1100.		0.0		[expr 2000+380.]
            theDomain.AddNode(node20601);

            var node20701 = new OpenSees.Components.NodeWrapper(20701, 6, 1100.0, 0.0, 2000 + 560.0);//node	20701		1100.		0.0		[expr 2000+560.]
            theDomain.AddNode(node20701);

            var node20801 = new OpenSees.Components.NodeWrapper(20801, 6, 1100.0, 0.0, 2000 + 740.0);//node	20801		1100.		0.0		[expr 2000+740.]
            theDomain.AddNode(node20801);

            var node20901 = new OpenSees.Components.NodeWrapper(20901, 6, 1100.0, 0.0, 2000 + 752.5);//node	20901		1100.		0.0		[expr 2000+752.5]
            theDomain.AddNode(node20901);

            var node21001 = new OpenSees.Components.NodeWrapper(21001, 6, 1100.0, 0.0, 2000 + 765.0);//node	21001		1100.		0.0		[expr 2000+765.]
            theDomain.AddNode(node21001);

            var node21101 = new OpenSees.Components.NodeWrapper(21101, 6, 1100.0, 0.0, 2000 + 777.5);//node	21101		1100.		0.0		[expr 2000+777.5]
            theDomain.AddNode(node21101);

            var node21201 = new OpenSees.Components.NodeWrapper(21201, 6, 1100.0, 0.0, 2000 + 790.0);//node	21201		1100.		0.0		[expr 2000+790.]
            theDomain.AddNode(node21201);

            var node21301 = new OpenSees.Components.NodeWrapper(21301, 6, 1100.0, 0.0, 2000 + 840.0);//node	21301		1100.		0.0		[expr 2000+840.]
            theDomain.AddNode(node21301);

            var node21401 = new OpenSees.Components.NodeWrapper(21401, 6, 1100.0, 0.0, 2000 + 1120.0);//node	21401		1100.		0.0		[expr 2000+1120.]
            theDomain.AddNode(node21401);

            var node21501 = new OpenSees.Components.NodeWrapper(21501, 6, 1100.0, 0.0, 2000 + 1130.0);//node	21501		1100.		0.0		[expr 2000+1130.]
            theDomain.AddNode(node21501);

            var node21601 = new OpenSees.Components.NodeWrapper(21601, 6, 1100.0, 0.0, 2000 + 1142.5);//node	21601		1100.		0.0		[expr 2000+1142.5]
            theDomain.AddNode(node21601);

            var node21701 = new OpenSees.Components.NodeWrapper(21701, 6, 1100.0, 0.0, 2000 + 1155.0);//node	21701		1100.		0.0		[expr 2000+1155.]
            theDomain.AddNode(node21701);

            var node21801 = new OpenSees.Components.NodeWrapper(21801, 6, 1100.0, 0.0, 2000 + 1167.5);//node	21801		1100.		0.0		[expr 2000+1167.5]
            theDomain.AddNode(node21801);

            var node21901 = new OpenSees.Components.NodeWrapper(21901, 6, 1100.0, 0.0, 2000 + 1180.0);//node	21901		1100.		0.0		[expr 2000+1180.]
            theDomain.AddNode(node21901);

            var node22001 = new OpenSees.Components.NodeWrapper(22001, 6, 1100.0, 0.0, 2000 + 1400.0);//node	22001		1100.		0.0		[expr 2000+1400.]
            theDomain.AddNode(node22001);

            var node30101 = new OpenSees.Components.NodeWrapper(30101, 6, 880.0, 0.0, 2000 + 0.0);//node	30101		880.		0.0		[expr 2000+0.]
            theDomain.AddNode(node30101);

            var node30201 = new OpenSees.Components.NodeWrapper(30201, 6, 863.75, 0.0, 2000 + 0.0);//node	30201		863.75		0.0		[expr 2000+0.]
            theDomain.AddNode(node30201);

            var node30301 = new OpenSees.Components.NodeWrapper(30301, 6, 852.5, 0.0, 2000 + 0.0);//node	30301		852.5		0.0		[expr 2000+0.]
            theDomain.AddNode(node30301);

            var node30401 = new OpenSees.Components.NodeWrapper(30401, 6, 851.25, 0.0, 2000 + 0.0);//node	30401		851.25		0.0		[expr 2000+0.]
            theDomain.AddNode(node30401);

            var node30501 = new OpenSees.Components.NodeWrapper(30501, 6, 838.75, 0.0, 2000 + 0.0);//node	30501		838.75		0.0		[expr 2000+0.]
            theDomain.AddNode(node30501);

            var node30601 = new OpenSees.Components.NodeWrapper(30601, 6, 826.25, 0.0, 2000 + 0.0);//node	30601		826.25		0.0		[expr 2000+0.]
            theDomain.AddNode(node30601);

            var node30701 = new OpenSees.Components.NodeWrapper(30701, 6, 813.75, 0.0, 2000 + 0.0);//node	30701		813.75		0.0		[expr 2000+0.]
            theDomain.AddNode(node30701);

            var node30801 = new OpenSees.Components.NodeWrapper(30801, 6, 660.0, 0.0, 2000 + 0.0);//node	30801		660.		0.0		[expr 2000+0.]
            theDomain.AddNode(node30801);

            var node30901 = new OpenSees.Components.NodeWrapper(30901, 6, 451.25, 0.0, 2000 + 0.0);//node	30901		451.25		0.0		[expr 2000+0.]
            theDomain.AddNode(node30901);

            var node31001 = new OpenSees.Components.NodeWrapper(31001, 6, 440.0, 0.0, 2000 + 0.0);//node	31001		440.		0.0		[expr 2000+0.]
            theDomain.AddNode(node31001);

            var node31101 = new OpenSees.Components.NodeWrapper(31101, 6, 438.75, 0.0, 2000 + 0.0);//node	31101		438.75		0.0		[expr 2000+0.]
            theDomain.AddNode(node31101);

            var node31201 = new OpenSees.Components.NodeWrapper(31201, 6, 426.25, 0.0, 2000 + 0.0);//node	31201		426.25		0.0		[expr 2000+0.]
            theDomain.AddNode(node31201);

            var node31301 = new OpenSees.Components.NodeWrapper(31301, 6, 413.75, 0.0, 2000 + 0.0);//node	31301		413.75		0.0		[expr 2000+0.]
            theDomain.AddNode(node31301);

            var node31401 = new OpenSees.Components.NodeWrapper(31401, 6, 401.25, 0.0, 2000 + 0.0);//node	31401		401.25		0.0		[expr 2000+0.]
            theDomain.AddNode(node31401);

            var node31501 = new OpenSees.Components.NodeWrapper(31501, 6, 220.0, 0.0, 2000 + 0.0);//node	31501		220.		0.0		[expr 2000+0.]
            theDomain.AddNode(node31501);

            var node40101 = new OpenSees.Components.NodeWrapper(40101, 6, 880.0, 0.0, 2000 + 1400.0);//node	40101		880.		0.0		[expr 2000+1400.]
            theDomain.AddNode(node40101);

            var node40201 = new OpenSees.Components.NodeWrapper(40201, 6, 697.394, 0.0, 2000 + 1400.0);//node	40201		697.394		0.0		[expr 2000+1400.]
            theDomain.AddNode(node40201);

            var node40301 = new OpenSees.Components.NodeWrapper(40301, 6, 684.894, 0.0, 2000 + 1400.0);//node	40301		684.894		0.0		[expr 2000+1400.]
            theDomain.AddNode(node40301);

            var node40401 = new OpenSees.Components.NodeWrapper(40401, 6, 672.394, 0.0, 2000 + 1400.0);//node	40401		672.394		0.0		[expr 2000+1400.]
            theDomain.AddNode(node40401);

            var node40501 = new OpenSees.Components.NodeWrapper(40501, 6, 660.0, 0.0, 2000 + 1400.0);//node	40501		660.		0.0		[expr 2000+1400.]
            theDomain.AddNode(node40501);

            var node40601 = new OpenSees.Components.NodeWrapper(40601, 6, 659.894, 0.0, 2000 + 1400.0);//node	40601		659.894		0.0		[expr 2000+1400.]
            theDomain.AddNode(node40601);

            var node40701 = new OpenSees.Components.NodeWrapper(40701, 6, 647.5, 0.0, 2000 + 1400.0);//node	40701		647.5		0.0		[expr 2000+1400.]
            theDomain.AddNode(node40701);

            var node40801 = new OpenSees.Components.NodeWrapper(40801, 6, 440.0, 0.0, 2000 + 1400.0);//node	40801		440.		0.0		[expr 2000+1400.]
            theDomain.AddNode(node40801);

            var node40901 = new OpenSees.Components.NodeWrapper(40901, 6, 280.207, 0.0, 2000 + 1400.0);//node	40901		280.207		0.0		[expr 2000+1400.]
            theDomain.AddNode(node40901);

            var node41001 = new OpenSees.Components.NodeWrapper(41001, 6, 268.957, 0.0, 2000 + 1400.0);//node	41001		268.957		0.0		[expr 2000+1400.]
            theDomain.AddNode(node41001);

            var node41101 = new OpenSees.Components.NodeWrapper(41101, 6, 267.707, 0.0, 2000 + 1400.0);//node	41101		267.707		0.0		[expr 2000+1400.]
            theDomain.AddNode(node41101);

            var node41201 = new OpenSees.Components.NodeWrapper(41201, 6, 255.207, 0.0, 2000 + 1400.0);//node	41201		255.207		0.0		[expr 2000+1400.]
            theDomain.AddNode(node41201);

            var node41301 = new OpenSees.Components.NodeWrapper(41301, 6, 242.707, 0.0, 2000 + 1400.0);//node	41301		242.707		0.0		[expr 2000+1400.]
            theDomain.AddNode(node41301);

            var node41401 = new OpenSees.Components.NodeWrapper(41401, 6, 230.207, 0.0, 2000 + 1400.0);//node	41401		230.207		0.0		[expr 2000+1400.]
            theDomain.AddNode(node41401);

            var node41501 = new OpenSees.Components.NodeWrapper(41501, 6, 220.0, 0.0, 2000 + 1400.0);//node	41501		220.		0.0		[expr 2000+1400.]
            theDomain.AddNode(node41501);

            foreach (var tag in new[] { 3020, 3040, 3050, 3060, 3070, 3090, 3110, 3120, 3130, 3140, 1180, 1150, 1140, 1110, 1120, 1090, 1080, 1010, 1020, 1050, 1040, 1170, 1030, 1100, 1160 })
                for (var dof = 0; dof < 6; dof++)
                    theDomain.AddSP_Constraint(new OpenSees.Components.Constraints.SP_ConstraintWrapper(tag, dof, 0.0, true));

            {
                var Ccr = new OpenSees.MatrixWrapper(6, 6);
                Ccr.Zero();
                Ccr[0, 0] = 1.0;
                Ccr[1, 1] = 1.0;
                Ccr[2, 2] = 1.0;
                Ccr[3, 3] = 1.0;
                Ccr[4, 4] = 1.0;
                Ccr[5, 5] = 1.0;

                var rcDOF = new OpenSees.IDWrapper(6);
                rcDOF[0] = 0;
                rcDOF[1] = 1;
                rcDOF[2] = 2;
                rcDOF[3] = 3;
                rcDOF[4] = 4;
                rcDOF[5] = 5;
                theDomain.AddMP_Constraint(new OpenSees.Components.Constraints.MP_ConstraintWrapper(4040, 10001, Ccr, rcDOF, rcDOF));
            }

            {
                var Ccr = new OpenSees.MatrixWrapper(6, 6);
                Ccr.Zero();
                Ccr[0, 0] = 1.0;
                Ccr[1, 1] = 1.0;
                Ccr[2, 2] = 1.0;
                Ccr[3, 3] = 1.0;
                Ccr[4, 4] = 1.0;
                Ccr[5, 5] = 1.0;

                var rcDOF = new OpenSees.IDWrapper(6);
                rcDOF[0] = 0;
                rcDOF[1] = 1;
                rcDOF[2] = 2;
                rcDOF[3] = 3;
                rcDOF[4] = 4;
                rcDOF[5] = 5;
                theDomain.AddMP_Constraint(new OpenSees.Components.Constraints.MP_ConstraintWrapper(4130, 20001, Ccr, rcDOF, rcDOF));

            }
            #endregion

            #region elements
            var xzplaneVec = new OpenSees.VectorWrapper(new double[] { 0, 1, 0 });
            var lineartag_TimoshenkoCrdT = new OpenSees.Elements.CrdTransfs.LinearCrdTransf3dWrapper(xzplaneVec);
            var lineartag_VerticalCrdT = new OpenSees.Elements.CrdTransfs.LinearCrdTransf3dWrapper(xzplaneVec);
            var lineartag_HorizontalCrdT = new OpenSees.Elements.CrdTransfs.LinearCrdTransf3dWrapper(xzplaneVec);

            //element elasticBeamColumn 	1		100		101	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_down		$lineartag_Vertical
            var ele1 = new OpenSees.Elements.ElasticBeam3dWrapper(1, Area_Perimeter_long_left, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left_mod_down, 100, 101, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele1);

            //element elasticBeamColumn 	2		101		102	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_down		$lineartag_Vertical
            var ele2 = new OpenSees.Elements.ElasticBeam3dWrapper(2, Area_Perimeter_long_left, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left_mod_down, 101, 102, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele2);

            //element elasticBeamColumn 	3		102		103	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_down		$lineartag_Vertical
            var ele3 = new OpenSees.Elements.ElasticBeam3dWrapper(3, Area_Perimeter_long_left, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left_mod_down, 102, 103, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele3);

            //element elasticBeamColumn 	4		103		104	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_down		$lineartag_Vertical
            var ele4 = new OpenSees.Elements.ElasticBeam3dWrapper(4, Area_Perimeter_long_left, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left_mod_down, 103, 104, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele4);

            //element elasticBeamColumn 	5		104		105	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_down		$lineartag_Vertical
            var ele5 = new OpenSees.Elements.ElasticBeam3dWrapper(5, Area_Perimeter_long_left, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left_mod_down, 104, 105, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele5);

            //element elasticBeamColumn 	6		105		106	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_down		$lineartag_Vertical
            var ele6 = new OpenSees.Elements.ElasticBeam3dWrapper(6, Area_Perimeter_long_left, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left_mod_down, 105, 106, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele6);

            //element elasticBeamColumn 	7		106		107	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_down		$lineartag_Vertical
            var ele7 = new OpenSees.Elements.ElasticBeam3dWrapper(7, Area_Perimeter_long_left, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left_mod_down, 106, 107, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele7);

            //element elasticBeamColumn 	8		107		108	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_mid			$lineartag_Vertical
            var ele8 = new OpenSees.Elements.ElasticBeam3dWrapper(8, Area_Perimeter_long_left, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left_mod_mid, 107, 108, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele8);

            //element elasticBeamColumn 	9		108		109	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_mid			$lineartag_Vertical
            var ele9 = new OpenSees.Elements.ElasticBeam3dWrapper(9, Area_Perimeter_long_left, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left_mod_mid, 108, 109, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele9);

            //element elasticBeamColumn 	129		109		110	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_mid			$lineartag_Vertical
            var ele129 = new OpenSees.Elements.ElasticBeam3dWrapper(129, Area_Perimeter_long_left, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left_mod_mid, 109, 110, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele129);

            //element elasticBeamColumn 	11		110		111	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_mid			$lineartag_Vertical
            var ele11 = new OpenSees.Elements.ElasticBeam3dWrapper(11, Area_Perimeter_long_left, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left_mod_mid, 110, 111, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele11);

            //element elasticBeamColumn 	12		111		112	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_mid			$lineartag_Vertical
            var ele12 = new OpenSees.Elements.ElasticBeam3dWrapper(12, Area_Perimeter_long_left, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left_mod_mid, 111, 112, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele12);

            //element elasticBeamColumn 	13		112		113	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_mid			$lineartag_Vertical
            var ele13 = new OpenSees.Elements.ElasticBeam3dWrapper(13, Area_Perimeter_long_left, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left_mod_mid, 112, 113, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele13);

            //element elasticBeamColumn 	14		113		114	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_up			$lineartag_Vertical
            var ele14 = new OpenSees.Elements.ElasticBeam3dWrapper(14, Area_Perimeter_long_left, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left_mod_up, 113, 114, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele14);

            //element elasticBeamColumn 	15		114		115	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_up			$lineartag_Vertical
            var ele15 = new OpenSees.Elements.ElasticBeam3dWrapper(15, Area_Perimeter_long_left, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left_mod_up, 114, 115, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele15);

            //element elasticBeamColumn 	16		115		116	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_up			$lineartag_Vertical
            var ele16 = new OpenSees.Elements.ElasticBeam3dWrapper(16, Area_Perimeter_long_left, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left_mod_up, 115, 116, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele16);

            //element elasticBeamColumn 	17		116		117	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_up			$lineartag_Vertical
            var ele17 = new OpenSees.Elements.ElasticBeam3dWrapper(17, Area_Perimeter_long_left, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left_mod_up, 116, 117, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele17);

            //element elasticBeamColumn 	18		117		118	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_up			$lineartag_Vertical
            var ele18 = new OpenSees.Elements.ElasticBeam3dWrapper(18, Area_Perimeter_long_left, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left_mod_up, 117, 118, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele18);

            //element elasticBeamColumn 	19		118		119	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_up			$lineartag_Vertical
            var ele19 = new OpenSees.Elements.ElasticBeam3dWrapper(19, Area_Perimeter_long_left, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left_mod_up, 118, 119, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele19);

            //element elasticBeamColumn 	130		119		120	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_down		$lineartag_Vertical
            var ele130 = new OpenSees.Elements.ElasticBeam3dWrapper(130, Area_Perimeter_long_left, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left_mod_down, 119, 120, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele130);

            //element elasticBeamColumn 	20		200		201	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_down		$lineartag_Vertical
            var ele20 = new OpenSees.Elements.ElasticBeam3dWrapper(20, Area_Perimeter_long_right, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right_mod_down, 200, 201, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele20);

            //element elasticBeamColumn 	21		201		202	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_down		$lineartag_Vertical
            var ele21 = new OpenSees.Elements.ElasticBeam3dWrapper(21, Area_Perimeter_long_right, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right_mod_down, 201, 202, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele21);

            //element elasticBeamColumn 	22		202		203	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_down		$lineartag_Vertical
            var ele22 = new OpenSees.Elements.ElasticBeam3dWrapper(22, Area_Perimeter_long_right, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right_mod_down, 202, 203, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele22);

            //element elasticBeamColumn 	23		203		204	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_down		$lineartag_Vertical
            var ele23 = new OpenSees.Elements.ElasticBeam3dWrapper(23, Area_Perimeter_long_right, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right_mod_down, 203, 204, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele23);

            //element elasticBeamColumn 	24		204		205	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_down		$lineartag_Vertical
            var ele24 = new OpenSees.Elements.ElasticBeam3dWrapper(24, Area_Perimeter_long_right, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right_mod_down, 204, 205, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele24);

            //element elasticBeamColumn 	25		205		206	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_down		$lineartag_Vertical
            var ele25 = new OpenSees.Elements.ElasticBeam3dWrapper(25, Area_Perimeter_long_right, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right_mod_down, 205, 206, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele25);

            //element elasticBeamColumn 	26		206		207	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_down		$lineartag_Vertical
            var ele26 = new OpenSees.Elements.ElasticBeam3dWrapper(26, Area_Perimeter_long_right, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right_mod_down, 206, 207, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele26);

            //element elasticBeamColumn 	27		207		208	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_mid		$lineartag_Vertical
            var ele27 = new OpenSees.Elements.ElasticBeam3dWrapper(27, Area_Perimeter_long_right, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right_mod_mid, 207, 208, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele27);

            //element elasticBeamColumn 	28		208		209	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_mid		$lineartag_Vertical
            var ele28 = new OpenSees.Elements.ElasticBeam3dWrapper(28, Area_Perimeter_long_right, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right_mod_mid, 208, 209, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele28);

            //element elasticBeamColumn 	29		209		210	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_mid		$lineartag_Vertical
            var ele29 = new OpenSees.Elements.ElasticBeam3dWrapper(29, Area_Perimeter_long_right, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right_mod_mid, 209, 210, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele29);

            //element elasticBeamColumn 	30		210		211	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_mid		$lineartag_Vertical
            var ele30 = new OpenSees.Elements.ElasticBeam3dWrapper(30, Area_Perimeter_long_right, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right_mod_mid, 210, 211, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele30);

            //element elasticBeamColumn 	31		211		212	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_mid		$lineartag_Vertical
            var ele31 = new OpenSees.Elements.ElasticBeam3dWrapper(31, Area_Perimeter_long_right, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right_mod_mid, 211, 212, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele31);

            //element elasticBeamColumn 	32		212		213	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_mid		$lineartag_Vertical
            var ele32 = new OpenSees.Elements.ElasticBeam3dWrapper(32, Area_Perimeter_long_right, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right_mod_mid, 212, 213, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele32);

            //element elasticBeamColumn 	33		213		214	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_mid		$lineartag_Vertical
            var ele33 = new OpenSees.Elements.ElasticBeam3dWrapper(33, Area_Perimeter_long_right, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right_mod_mid, 213, 214, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele33);

            //element elasticBeamColumn 	34		214		215	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_up			$lineartag_Vertical
            var ele34 = new OpenSees.Elements.ElasticBeam3dWrapper(34, Area_Perimeter_long_right, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right_mod_up, 214, 215, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele34);

            //element elasticBeamColumn 	35		215		216	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_up			$lineartag_Vertical
            var ele35 = new OpenSees.Elements.ElasticBeam3dWrapper(35, Area_Perimeter_long_right, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right_mod_up, 215, 216, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele35);

            //element elasticBeamColumn 	36		216		217	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_up			$lineartag_Vertical
            var ele36 = new OpenSees.Elements.ElasticBeam3dWrapper(36, Area_Perimeter_long_right, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right_mod_up, 216, 217, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele36);

            //element elasticBeamColumn 	37		217		218	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_up			$lineartag_Vertical
            var ele37 = new OpenSees.Elements.ElasticBeam3dWrapper(37, Area_Perimeter_long_right, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right_mod_up, 217, 218, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele37);

            //element elasticBeamColumn 	38		218		219	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_up			$lineartag_Vertical
            var ele38 = new OpenSees.Elements.ElasticBeam3dWrapper(38, Area_Perimeter_long_right, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right_mod_up, 218, 219, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele38);

            //element elasticBeamColumn 	131		219		220	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_up			$lineartag_Vertical
            var ele131 = new OpenSees.Elements.ElasticBeam3dWrapper(131, Area_Perimeter_long_right, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right_mod_up, 219, 220, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele131);

            //element elasticBeamColumn 	200		10001		10101	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_down		$lineartag_Vertical
            var ele200 = new OpenSees.Elements.ElasticBeam3dWrapper(200, Area_Perimeter_long_left, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left_mod_down, 10001, 10101, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele200);

            //element elasticBeamColumn 	201		10101		10201	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_down		$lineartag_Vertical
            var ele201 = new OpenSees.Elements.ElasticBeam3dWrapper(201, Area_Perimeter_long_left, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left_mod_down, 10101, 10201, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele201);

            //element elasticBeamColumn 	202		10201		10301	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_down		$lineartag_Vertical
            var ele202 = new OpenSees.Elements.ElasticBeam3dWrapper(202, Area_Perimeter_long_left, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left_mod_down, 10201, 10301, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele202);

            //element elasticBeamColumn 	203		10301		10401	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_down		$lineartag_Vertical
            var ele203 = new OpenSees.Elements.ElasticBeam3dWrapper(203, Area_Perimeter_long_left, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left_mod_down, 10301, 10401, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele203);

            //element elasticBeamColumn 	204		10401		10501	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_down		$lineartag_Vertical
            var ele204 = new OpenSees.Elements.ElasticBeam3dWrapper(204, Area_Perimeter_long_left, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left_mod_down, 10401, 10501, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele204);

            //element elasticBeamColumn 	205		10501		10601	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_down		$lineartag_Vertical
            var ele205 = new OpenSees.Elements.ElasticBeam3dWrapper(205, Area_Perimeter_long_left, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left_mod_down, 10501, 10601, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele205);

            //element elasticBeamColumn 	206		10601		10701	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_down		$lineartag_Vertical
            var ele206 = new OpenSees.Elements.ElasticBeam3dWrapper(206, Area_Perimeter_long_left, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left_mod_down, 10601, 10701, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele206);

            //element elasticBeamColumn 	207		10701		10801	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_mid			$lineartag_Vertical
            var ele207 = new OpenSees.Elements.ElasticBeam3dWrapper(207, Area_Perimeter_long_left, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left_mod_mid, 10701, 10801, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele207);

            //element elasticBeamColumn 	208		10801		10901	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_mid			$lineartag_Vertical
            var ele208 = new OpenSees.Elements.ElasticBeam3dWrapper(208, Area_Perimeter_long_left, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left_mod_mid, 10801, 10901, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele208);

            //element elasticBeamColumn 	209		10901		11001	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_mid			$lineartag_Vertical
            var ele209 = new OpenSees.Elements.ElasticBeam3dWrapper(209, Area_Perimeter_long_left, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left_mod_mid, 10901, 11001, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele209);

            //element elasticBeamColumn 	210		11001		11101	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_mid			$lineartag_Vertical
            var ele210 = new OpenSees.Elements.ElasticBeam3dWrapper(210, Area_Perimeter_long_left, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left_mod_mid, 11001, 11101, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele210);

            //element elasticBeamColumn 	211		11101		11201	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_mid			$lineartag_Vertical
            var ele211 = new OpenSees.Elements.ElasticBeam3dWrapper(211, Area_Perimeter_long_left, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left_mod_mid, 11101, 11201, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele211);

            //element elasticBeamColumn 	212		11201		11301	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_mid			$lineartag_Vertical
            var ele212 = new OpenSees.Elements.ElasticBeam3dWrapper(212, Area_Perimeter_long_left, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left_mod_mid, 11201, 11301, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele212);

            //element elasticBeamColumn 	213		11301		11401	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_up			$lineartag_Vertical
            var ele213 = new OpenSees.Elements.ElasticBeam3dWrapper(213, Area_Perimeter_long_left, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left_mod_up, 11301, 11401, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele213);

            //element elasticBeamColumn 	214		11401		11501	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_up			$lineartag_Vertical
            var ele214 = new OpenSees.Elements.ElasticBeam3dWrapper(214, Area_Perimeter_long_left, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left_mod_up, 11401, 11501, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele214);

            //element elasticBeamColumn 	215		11501		11601	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_up			$lineartag_Vertical
            var ele215 = new OpenSees.Elements.ElasticBeam3dWrapper(215, Area_Perimeter_long_left, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left_mod_up, 11501, 11601, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele215);

            //element elasticBeamColumn 	216		11601		11701	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_up			$lineartag_Vertical
            var ele216 = new OpenSees.Elements.ElasticBeam3dWrapper(216, Area_Perimeter_long_left, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left_mod_up, 11601, 11701, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele216);

            //element elasticBeamColumn 	217		11701		11801	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_up			$lineartag_Vertical
            var ele217 = new OpenSees.Elements.ElasticBeam3dWrapper(217, Area_Perimeter_long_left, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left_mod_up, 11701, 11801, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele217);

            //element elasticBeamColumn 	218		11801		11901	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_up			$lineartag_Vertical
            var ele218 = new OpenSees.Elements.ElasticBeam3dWrapper(218, Area_Perimeter_long_left, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left_mod_up, 11801, 11901, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele218);

            //element elasticBeamColumn 	219		11901		12001	 		$Area_Perimeter_long_left 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_left	$Iy_Perimeter_long_left			$Iz_Perimeter_long_left_mod_down		$lineartag_Vertical
            var ele219 = new OpenSees.Elements.ElasticBeam3dWrapper(219, Area_Perimeter_long_left, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_left, Iy_Perimeter_long_left, Iz_Perimeter_long_left_mod_down, 11901, 12001, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele219);

            //element elasticBeamColumn 	220		20001		20101	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_down		$lineartag_Vertical
            var ele220 = new OpenSees.Elements.ElasticBeam3dWrapper(220, Area_Perimeter_long_right, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right_mod_down, 20001, 20101, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele220);

            //element elasticBeamColumn 	221		20101		20201	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_down		$lineartag_Vertical
            var ele221 = new OpenSees.Elements.ElasticBeam3dWrapper(221, Area_Perimeter_long_right, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right_mod_down, 20101, 20201, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele221);

            //element elasticBeamColumn 	222		20201		20301	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_down		$lineartag_Vertical
            var ele222 = new OpenSees.Elements.ElasticBeam3dWrapper(222, Area_Perimeter_long_right, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right_mod_down, 20201, 20301, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele222);

            //element elasticBeamColumn 	223		20301		20401	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_down		$lineartag_Vertical
            var ele223 = new OpenSees.Elements.ElasticBeam3dWrapper(223, Area_Perimeter_long_right, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right_mod_down, 20301, 20401, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele223);

            //element elasticBeamColumn 	224		20401		20501	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_down		$lineartag_Vertical
            var ele224 = new OpenSees.Elements.ElasticBeam3dWrapper(224, Area_Perimeter_long_right, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right_mod_down, 20401, 20501, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele224);

            //element elasticBeamColumn 	225		20501		20601	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_down		$lineartag_Vertical
            var ele225 = new OpenSees.Elements.ElasticBeam3dWrapper(225, Area_Perimeter_long_right, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right_mod_down, 20501, 20601, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele225);

            //element elasticBeamColumn 	226		20601		20701	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_down		$lineartag_Vertical
            var ele226 = new OpenSees.Elements.ElasticBeam3dWrapper(226, Area_Perimeter_long_right, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right_mod_down, 20601, 20701, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele226);

            //element elasticBeamColumn 	227		20701		20801	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_mid		$lineartag_Vertical
            var ele227 = new OpenSees.Elements.ElasticBeam3dWrapper(227, Area_Perimeter_long_right, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right_mod_mid, 20701, 20801, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele227);

            //element elasticBeamColumn 	228		20801		20901	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_mid		$lineartag_Vertical
            var ele228 = new OpenSees.Elements.ElasticBeam3dWrapper(228, Area_Perimeter_long_right, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right_mod_mid, 20801, 20901, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele228);

            //element elasticBeamColumn 	229		20901		21001	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_mid		$lineartag_Vertical
            var ele229 = new OpenSees.Elements.ElasticBeam3dWrapper(229, Area_Perimeter_long_right, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right_mod_mid, 20901, 21001, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele229);

            //element elasticBeamColumn 	230		21001		21101	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_mid		$lineartag_Vertical
            var ele230 = new OpenSees.Elements.ElasticBeam3dWrapper(230, Area_Perimeter_long_right, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right_mod_mid, 21001, 21101, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele230);

            //element elasticBeamColumn 	231		21101		21201	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_mid		$lineartag_Vertical
            var ele231 = new OpenSees.Elements.ElasticBeam3dWrapper(231, Area_Perimeter_long_right, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right_mod_mid, 21101, 21201, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele231);

            //element elasticBeamColumn 	232		21201		21301	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_mid		$lineartag_Vertical
            var ele232 = new OpenSees.Elements.ElasticBeam3dWrapper(232, Area_Perimeter_long_right, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right_mod_mid, 21201, 21301, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele232);

            //element elasticBeamColumn 	233		21301		21401	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_mid		$lineartag_Vertical
            var ele233 = new OpenSees.Elements.ElasticBeam3dWrapper(233, Area_Perimeter_long_right, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right_mod_mid, 21301, 21401, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele233);

            //element elasticBeamColumn 	234		21401		21501	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_up			$lineartag_Vertical
            var ele234 = new OpenSees.Elements.ElasticBeam3dWrapper(234, Area_Perimeter_long_right, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right_mod_up, 21401, 21501, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele234);

            //element elasticBeamColumn 	235		21501		21601	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_up			$lineartag_Vertical
            var ele235 = new OpenSees.Elements.ElasticBeam3dWrapper(235, Area_Perimeter_long_right, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right_mod_up, 21501, 21601, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele235);

            //element elasticBeamColumn 	236		21601		21701	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_up			$lineartag_Vertical
            var ele236 = new OpenSees.Elements.ElasticBeam3dWrapper(236, Area_Perimeter_long_right, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right_mod_up, 21601, 21701, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele236);

            //element elasticBeamColumn 	237		21701		21801	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_up			$lineartag_Vertical
            var ele237 = new OpenSees.Elements.ElasticBeam3dWrapper(237, Area_Perimeter_long_right, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right_mod_up, 21701, 21801, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele237);

            //element elasticBeamColumn 	238		21801		21901	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_up			$lineartag_Vertical
            var ele238 = new OpenSees.Elements.ElasticBeam3dWrapper(238, Area_Perimeter_long_right, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right_mod_up, 21801, 21901, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele238);

            //element elasticBeamColumn 	239		21901		22001	 		$Area_Perimeter_long_right 	$E_Timber_Long	$G_Timber_Long	$Torsional_J_Perimeter_long_right	$Iy_Perimeter_long_right		$Iz_Perimeter_long_right_mod_up			$lineartag_Vertical
            var ele239 = new OpenSees.Elements.ElasticBeam3dWrapper(239, Area_Perimeter_long_right, E_Timber_Long, G_Timber_Long, Torsional_J_Perimeter_long_right, Iy_Perimeter_long_right, Iz_Perimeter_long_right_mod_up, 21901, 22001, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele239);

            //element elasticBeamColumn 	39		200		301			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
            var ele39 = new OpenSees.Elements.ElasticBeam3dWrapper(39, Area_Perimeter_trans_down, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_down, Iy_Perimeter_trans_down, Iz_Perimeter_trans_down_mod, 200, 301, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele39);

            //element elasticBeamColumn 	40		301		302			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
            var ele40 = new OpenSees.Elements.ElasticBeam3dWrapper(40, Area_Perimeter_trans_down, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_down, Iy_Perimeter_trans_down, Iz_Perimeter_trans_down_mod, 301, 302, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele40);

            //element elasticBeamColumn 	41		302		303			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
            var ele41 = new OpenSees.Elements.ElasticBeam3dWrapper(41, Area_Perimeter_trans_down, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_down, Iy_Perimeter_trans_down, Iz_Perimeter_trans_down_mod, 302, 303, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele41);

            //element elasticBeamColumn 	42		303		304			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
            var ele42 = new OpenSees.Elements.ElasticBeam3dWrapper(42, Area_Perimeter_trans_down, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_down, Iy_Perimeter_trans_down, Iz_Perimeter_trans_down_mod, 303, 304, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele42);

            //element elasticBeamColumn 	43		304		305			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
            var ele43 = new OpenSees.Elements.ElasticBeam3dWrapper(43, Area_Perimeter_trans_down, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_down, Iy_Perimeter_trans_down, Iz_Perimeter_trans_down_mod, 304, 305, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele43);

            //element elasticBeamColumn 	44		305		306			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
            var ele44 = new OpenSees.Elements.ElasticBeam3dWrapper(44, Area_Perimeter_trans_down, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_down, Iy_Perimeter_trans_down, Iz_Perimeter_trans_down_mod, 305, 306, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele44);

            //element elasticBeamColumn 	45		306		307			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
            var ele45 = new OpenSees.Elements.ElasticBeam3dWrapper(45, Area_Perimeter_trans_down, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_down, Iy_Perimeter_trans_down, Iz_Perimeter_trans_down_mod, 306, 307, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele45);

            //element elasticBeamColumn 	46		307		308			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
            var ele46 = new OpenSees.Elements.ElasticBeam3dWrapper(46, Area_Perimeter_trans_down, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_down, Iy_Perimeter_trans_down, Iz_Perimeter_trans_down_mod, 307, 308, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele46);

            //element elasticBeamColumn 	47		308		309			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
            var ele47 = new OpenSees.Elements.ElasticBeam3dWrapper(47, Area_Perimeter_trans_down, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_down, Iy_Perimeter_trans_down, Iz_Perimeter_trans_down_mod, 308, 309, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele47);

            //element elasticBeamColumn 	48		309		310			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
            var ele48 = new OpenSees.Elements.ElasticBeam3dWrapper(48, Area_Perimeter_trans_down, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_down, Iy_Perimeter_trans_down, Iz_Perimeter_trans_down_mod, 309, 310, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele48);

            //element elasticBeamColumn 	49		310		311			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
            var ele49 = new OpenSees.Elements.ElasticBeam3dWrapper(49, Area_Perimeter_trans_down, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_down, Iy_Perimeter_trans_down, Iz_Perimeter_trans_down_mod, 310, 311, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele49);

            //element elasticBeamColumn 	50		311		312			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
            var ele50 = new OpenSees.Elements.ElasticBeam3dWrapper(50, Area_Perimeter_trans_down, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_down, Iy_Perimeter_trans_down, Iz_Perimeter_trans_down_mod, 311, 312, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele50);

            //element elasticBeamColumn 	51		312		313			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
            var ele51 = new OpenSees.Elements.ElasticBeam3dWrapper(51, Area_Perimeter_trans_down, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_down, Iy_Perimeter_trans_down, Iz_Perimeter_trans_down_mod, 312, 313, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele51);

            //element elasticBeamColumn 	52		313		314			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
            var ele52 = new OpenSees.Elements.ElasticBeam3dWrapper(52, Area_Perimeter_trans_down, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_down, Iy_Perimeter_trans_down, Iz_Perimeter_trans_down_mod, 313, 314, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele52);

            //element elasticBeamColumn 	53		314		315			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
            var ele53 = new OpenSees.Elements.ElasticBeam3dWrapper(53, Area_Perimeter_trans_down, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_down, Iy_Perimeter_trans_down, Iz_Perimeter_trans_down_mod, 314, 315, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele53);

            //element elasticBeamColumn 	54		315		100			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
            var ele54 = new OpenSees.Elements.ElasticBeam3dWrapper(54, Area_Perimeter_trans_down, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_down, Iy_Perimeter_trans_down, Iz_Perimeter_trans_down_mod, 315, 100, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele54);

            //element elasticBeamColumn 	55		220		401			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
            var ele55 = new OpenSees.Elements.ElasticBeam3dWrapper(55, Area_Perimeter_trans_up, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_up, Iy_Perimeter_trans_up, Iz_Perimeter_trans_up_mod, 220, 401, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele55);

            //element elasticBeamColumn 	56		401		402			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
            var ele56 = new OpenSees.Elements.ElasticBeam3dWrapper(56, Area_Perimeter_trans_up, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_up, Iy_Perimeter_trans_up, Iz_Perimeter_trans_up_mod, 401, 402, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele56);

            //element elasticBeamColumn 	57		402		403			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
            var ele57 = new OpenSees.Elements.ElasticBeam3dWrapper(57, Area_Perimeter_trans_up, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_up, Iy_Perimeter_trans_up, Iz_Perimeter_trans_up_mod, 402, 403, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele57);

            //element elasticBeamColumn 	58		403		404			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
            var ele58 = new OpenSees.Elements.ElasticBeam3dWrapper(58, Area_Perimeter_trans_up, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_up, Iy_Perimeter_trans_up, Iz_Perimeter_trans_up_mod, 403, 404, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele58);

            //element elasticBeamColumn 	59		404		405			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
            var ele59 = new OpenSees.Elements.ElasticBeam3dWrapper(59, Area_Perimeter_trans_up, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_up, Iy_Perimeter_trans_up, Iz_Perimeter_trans_up_mod, 404, 405, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele59);

            //element elasticBeamColumn 	60		405		406			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
            var ele60 = new OpenSees.Elements.ElasticBeam3dWrapper(60, Area_Perimeter_trans_up, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_up, Iy_Perimeter_trans_up, Iz_Perimeter_trans_up_mod, 405, 406, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele60);

            //element elasticBeamColumn 	61		406		407			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
            var ele61 = new OpenSees.Elements.ElasticBeam3dWrapper(61, Area_Perimeter_trans_up, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_up, Iy_Perimeter_trans_up, Iz_Perimeter_trans_up_mod, 406, 407, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele61);

            //element elasticBeamColumn 	62		407		408			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
            var ele62 = new OpenSees.Elements.ElasticBeam3dWrapper(62, Area_Perimeter_trans_up, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_up, Iy_Perimeter_trans_up, Iz_Perimeter_trans_up_mod, 407, 408, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele62);

            //element elasticBeamColumn 	63		408		409			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
            var ele63 = new OpenSees.Elements.ElasticBeam3dWrapper(63, Area_Perimeter_trans_up, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_up, Iy_Perimeter_trans_up, Iz_Perimeter_trans_up_mod, 408, 409, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele63);

            //element elasticBeamColumn 	64		409		410			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
            var ele64 = new OpenSees.Elements.ElasticBeam3dWrapper(64, Area_Perimeter_trans_up, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_up, Iy_Perimeter_trans_up, Iz_Perimeter_trans_up_mod, 409, 410, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele64);

            //element elasticBeamColumn 	65		410		411			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
            var ele65 = new OpenSees.Elements.ElasticBeam3dWrapper(65, Area_Perimeter_trans_up, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_up, Iy_Perimeter_trans_up, Iz_Perimeter_trans_up_mod, 410, 411, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele65);

            //element elasticBeamColumn 	66		411		412			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
            var ele66 = new OpenSees.Elements.ElasticBeam3dWrapper(66, Area_Perimeter_trans_up, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_up, Iy_Perimeter_trans_up, Iz_Perimeter_trans_up_mod, 411, 412, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele66);

            //element elasticBeamColumn 	67		412		413			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
            var ele67 = new OpenSees.Elements.ElasticBeam3dWrapper(67, Area_Perimeter_trans_up, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_up, Iy_Perimeter_trans_up, Iz_Perimeter_trans_up_mod, 412, 413, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele67);

            //element elasticBeamColumn 	68		413		414			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
            var ele68 = new OpenSees.Elements.ElasticBeam3dWrapper(68, Area_Perimeter_trans_up, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_up, Iy_Perimeter_trans_up, Iz_Perimeter_trans_up_mod, 413, 414, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele68);

            //element elasticBeamColumn 	69		414		415			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
            var ele69 = new OpenSees.Elements.ElasticBeam3dWrapper(69, Area_Perimeter_trans_up, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_up, Iy_Perimeter_trans_up, Iz_Perimeter_trans_up_mod, 414, 415, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele69);

            //element elasticBeamColumn 	70		415		120			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
            var ele70 = new OpenSees.Elements.ElasticBeam3dWrapper(70, Area_Perimeter_trans_up, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_up, Iy_Perimeter_trans_up, Iz_Perimeter_trans_up_mod, 415, 120, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele70);

            //element elasticBeamColumn 	250		20001		30101			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
            var ele250 = new OpenSees.Elements.ElasticBeam3dWrapper(250, Area_Perimeter_trans_down, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_down, Iy_Perimeter_trans_down, Iz_Perimeter_trans_down_mod, 20001, 30101, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele250);

            //element elasticBeamColumn 	251		30101		30201			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
            var ele251 = new OpenSees.Elements.ElasticBeam3dWrapper(251, Area_Perimeter_trans_down, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_down, Iy_Perimeter_trans_down, Iz_Perimeter_trans_down_mod, 30101, 30201, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele251);

            //element elasticBeamColumn 	252		30201		30301			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
            var ele252 = new OpenSees.Elements.ElasticBeam3dWrapper(252, Area_Perimeter_trans_down, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_down, Iy_Perimeter_trans_down, Iz_Perimeter_trans_down_mod, 30201, 30301, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele252);

            //element elasticBeamColumn 	253		30301		30401			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
            var ele253 = new OpenSees.Elements.ElasticBeam3dWrapper(253, Area_Perimeter_trans_down, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_down, Iy_Perimeter_trans_down, Iz_Perimeter_trans_down_mod, 30301, 30401, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele253);

            //element elasticBeamColumn 	254		30401		30501			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
            var ele254 = new OpenSees.Elements.ElasticBeam3dWrapper(254, Area_Perimeter_trans_down, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_down, Iy_Perimeter_trans_down, Iz_Perimeter_trans_down_mod, 30401, 30501, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele254);

            //element elasticBeamColumn 	255		30501		30601			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
            var ele255 = new OpenSees.Elements.ElasticBeam3dWrapper(255, Area_Perimeter_trans_down, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_down, Iy_Perimeter_trans_down, Iz_Perimeter_trans_down_mod, 30501, 30601, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele255);

            //element elasticBeamColumn 	256		30601		30701			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
            var ele256 = new OpenSees.Elements.ElasticBeam3dWrapper(256, Area_Perimeter_trans_down, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_down, Iy_Perimeter_trans_down, Iz_Perimeter_trans_down_mod, 30601, 30701, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele256);

            //element elasticBeamColumn 	257		30701		30801			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
            var ele257 = new OpenSees.Elements.ElasticBeam3dWrapper(257, Area_Perimeter_trans_down, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_down, Iy_Perimeter_trans_down, Iz_Perimeter_trans_down_mod, 30701, 30801, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele257);

            //element elasticBeamColumn 	258		30801		30901			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
            var ele258 = new OpenSees.Elements.ElasticBeam3dWrapper(258, Area_Perimeter_trans_down, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_down, Iy_Perimeter_trans_down, Iz_Perimeter_trans_down_mod, 30801, 30901, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele258);

            //element elasticBeamColumn 	259		30901		31001			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
            var ele259 = new OpenSees.Elements.ElasticBeam3dWrapper(259, Area_Perimeter_trans_down, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_down, Iy_Perimeter_trans_down, Iz_Perimeter_trans_down_mod, 30901, 31001, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele259);

            //element elasticBeamColumn 	260		31001		31101			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
            var ele260 = new OpenSees.Elements.ElasticBeam3dWrapper(260, Area_Perimeter_trans_down, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_down, Iy_Perimeter_trans_down, Iz_Perimeter_trans_down_mod, 31001, 31101, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele260);

            //element elasticBeamColumn 	261		31101		31201			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
            var ele261 = new OpenSees.Elements.ElasticBeam3dWrapper(261, Area_Perimeter_trans_down, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_down, Iy_Perimeter_trans_down, Iz_Perimeter_trans_down_mod, 31101, 31201, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele261);

            //element elasticBeamColumn 	262		31201		31301			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
            var ele262 = new OpenSees.Elements.ElasticBeam3dWrapper(262, Area_Perimeter_trans_down, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_down, Iy_Perimeter_trans_down, Iz_Perimeter_trans_down_mod, 31201, 31301, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele262);

            //element elasticBeamColumn 	263		31301		31401			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
            var ele263 = new OpenSees.Elements.ElasticBeam3dWrapper(263, Area_Perimeter_trans_down, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_down, Iy_Perimeter_trans_down, Iz_Perimeter_trans_down_mod, 31301, 31401, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele263);

            //element elasticBeamColumn 	264		31401		31501			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
            var ele264 = new OpenSees.Elements.ElasticBeam3dWrapper(264, Area_Perimeter_trans_down, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_down, Iy_Perimeter_trans_down, Iz_Perimeter_trans_down_mod, 31401, 31501, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele264);

            //element elasticBeamColumn 	265		31501		10001			$Area_Perimeter_trans_down	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_down 	$Iy_Perimeter_trans_down	 $Iz_Perimeter_trans_down_mod	 $lineartag_Horizontal
            var ele265 = new OpenSees.Elements.ElasticBeam3dWrapper(265, Area_Perimeter_trans_down, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_down, Iy_Perimeter_trans_down, Iz_Perimeter_trans_down_mod, 31501, 10001, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele265);

            //element elasticBeamColumn 	266		22001		40101			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
            var ele266 = new OpenSees.Elements.ElasticBeam3dWrapper(266, Area_Perimeter_trans_up, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_up, Iy_Perimeter_trans_up, Iz_Perimeter_trans_up_mod, 22001, 40101, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele266);

            //element elasticBeamColumn 	267		40101		40201			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
            var ele267 = new OpenSees.Elements.ElasticBeam3dWrapper(267, Area_Perimeter_trans_up, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_up, Iy_Perimeter_trans_up, Iz_Perimeter_trans_up_mod, 40101, 40201, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele267);

            //element elasticBeamColumn 	268		40201		40301			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
            var ele268 = new OpenSees.Elements.ElasticBeam3dWrapper(268, Area_Perimeter_trans_up, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_up, Iy_Perimeter_trans_up, Iz_Perimeter_trans_up_mod, 40201, 40301, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele268);

            //element elasticBeamColumn 	269		40301		40401			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
            var ele269 = new OpenSees.Elements.ElasticBeam3dWrapper(269, Area_Perimeter_trans_up, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_up, Iy_Perimeter_trans_up, Iz_Perimeter_trans_up_mod, 40301, 40401, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele269);

            //element elasticBeamColumn 	270		40401		40501			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
            var ele270 = new OpenSees.Elements.ElasticBeam3dWrapper(270, Area_Perimeter_trans_up, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_up, Iy_Perimeter_trans_up, Iz_Perimeter_trans_up_mod, 40401, 40501, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele270);

            //element elasticBeamColumn 	271		40501		40601			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
            var ele271 = new OpenSees.Elements.ElasticBeam3dWrapper(271, Area_Perimeter_trans_up, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_up, Iy_Perimeter_trans_up, Iz_Perimeter_trans_up_mod, 40501, 40601, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele271);

            //element elasticBeamColumn 	272		40601		40701			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
            var ele272 = new OpenSees.Elements.ElasticBeam3dWrapper(272, Area_Perimeter_trans_up, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_up, Iy_Perimeter_trans_up, Iz_Perimeter_trans_up_mod, 40601, 40701, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele272);

            //element elasticBeamColumn 	273		40701		40801			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
            var ele273 = new OpenSees.Elements.ElasticBeam3dWrapper(273, Area_Perimeter_trans_up, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_up, Iy_Perimeter_trans_up, Iz_Perimeter_trans_up_mod, 40701, 40801, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele273);

            //element elasticBeamColumn 	274		40801		40901			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
            var ele274 = new OpenSees.Elements.ElasticBeam3dWrapper(274, Area_Perimeter_trans_up, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_up, Iy_Perimeter_trans_up, Iz_Perimeter_trans_up_mod, 40801, 40901, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele274);

            //element elasticBeamColumn 	275		40901		41001			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
            var ele275 = new OpenSees.Elements.ElasticBeam3dWrapper(275, Area_Perimeter_trans_up, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_up, Iy_Perimeter_trans_up, Iz_Perimeter_trans_up_mod, 40901, 41001, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele275);

            //element elasticBeamColumn 	276		41001		41101			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
            var ele276 = new OpenSees.Elements.ElasticBeam3dWrapper(276, Area_Perimeter_trans_up, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_up, Iy_Perimeter_trans_up, Iz_Perimeter_trans_up_mod, 41001, 41101, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele276);

            //element elasticBeamColumn 	277		41101		41201			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
            var ele277 = new OpenSees.Elements.ElasticBeam3dWrapper(277, Area_Perimeter_trans_up, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_up, Iy_Perimeter_trans_up, Iz_Perimeter_trans_up_mod, 41101, 41201, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele277);

            //element elasticBeamColumn 	278		41201		41301			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
            var ele278 = new OpenSees.Elements.ElasticBeam3dWrapper(278, Area_Perimeter_trans_up, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_up, Iy_Perimeter_trans_up, Iz_Perimeter_trans_up_mod, 41201, 41301, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele278);

            //element elasticBeamColumn 	279		41301		41401			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
            var ele279 = new OpenSees.Elements.ElasticBeam3dWrapper(279, Area_Perimeter_trans_up, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_up, Iy_Perimeter_trans_up, Iz_Perimeter_trans_up_mod, 41301, 41401, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele279);

            //element elasticBeamColumn 	280		41401		41501			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
            var ele280 = new OpenSees.Elements.ElasticBeam3dWrapper(280, Area_Perimeter_trans_up, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_up, Iy_Perimeter_trans_up, Iz_Perimeter_trans_up_mod, 41401, 41501, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele280);

            //element elasticBeamColumn 	281		41501		12001			$Area_Perimeter_trans_up	$E_Timber_Trans		$G_Timber_Trans	 $Torsional_J_Perimeter_trans_up 	$Iy_Perimeter_trans_up	 	$Iz_Perimeter_trans_up_mod		$lineartag_Horizontal
            var ele281 = new OpenSees.Elements.ElasticBeam3dWrapper(281, Area_Perimeter_trans_up, E_Timber_Trans, G_Timber_Trans, Torsional_J_Perimeter_trans_up, Iy_Perimeter_trans_up, Iz_Perimeter_trans_up_mod, 41501, 12001, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele281);

            //element elasticBeamColumn 71	315		415		 $Area_inside_long		$E_Timber_Long		$G_Timber_Long		$Torsional_J_inside_long		$Iy_inside_long		$Iz_inside_long_mod		$lineartag_Vertical
            var ele71 = new OpenSees.Elements.ElasticBeam3dWrapper(71, Area_inside_long, E_Timber_Long, G_Timber_Long, Torsional_J_inside_long, Iy_inside_long, Iz_inside_long_mod, 315, 415, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele71);

            //element elasticBeamColumn 72	310		408		 $Area_inside_long		$E_Timber_Long		$G_Timber_Long		$Torsional_J_inside_long		$Iy_inside_long		$Iz_inside_long_mod		$lineartag_Vertical
            var ele72 = new OpenSees.Elements.ElasticBeam3dWrapper(72, Area_inside_long, E_Timber_Long, G_Timber_Long, Torsional_J_inside_long, Iy_inside_long, Iz_inside_long_mod, 310, 408, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele72);

            //element elasticBeamColumn 73	308		405		 $Area_inside_long		$E_Timber_Long		$G_Timber_Long		$Torsional_J_inside_long		$Iy_inside_long		$Iz_inside_long_mod		$lineartag_Vertical
            var ele73 = new OpenSees.Elements.ElasticBeam3dWrapper(73, Area_inside_long, E_Timber_Long, G_Timber_Long, Torsional_J_inside_long, Iy_inside_long, Iz_inside_long_mod, 308, 405, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele73);

            //element elasticBeamColumn 74	301		401		 $Area_inside_long		$E_Timber_Long		$G_Timber_Long		$Torsional_J_inside_long		$Iy_inside_long		$Iz_inside_long_mod		$lineartag_Vertical
            var ele74 = new OpenSees.Elements.ElasticBeam3dWrapper(74, Area_inside_long, E_Timber_Long, G_Timber_Long, Torsional_J_inside_long, Iy_inside_long, Iz_inside_long_mod, 301, 401, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele74);

            //element elasticBeamColumn 300	31501		41501		 $Area_inside_long		$E_Timber_Long		$G_Timber_Long		$Torsional_J_inside_long		$Iy_inside_long		$Iz_inside_long_mod		$lineartag_Vertical
            var ele300 = new OpenSees.Elements.ElasticBeam3dWrapper(300, Area_inside_long, E_Timber_Long, G_Timber_Long, Torsional_J_inside_long, Iy_inside_long, Iz_inside_long_mod, 31501, 41501, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele300);

            //element elasticBeamColumn 301	31001		40801		 $Area_inside_long		$E_Timber_Long		$G_Timber_Long		$Torsional_J_inside_long		$Iy_inside_long		$Iz_inside_long_mod		$lineartag_Vertical
            var ele301 = new OpenSees.Elements.ElasticBeam3dWrapper(301, Area_inside_long, E_Timber_Long, G_Timber_Long, Torsional_J_inside_long, Iy_inside_long, Iz_inside_long_mod, 31001, 40801, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele301);

            //element elasticBeamColumn 302	30801		40501		 $Area_inside_long		$E_Timber_Long		$G_Timber_Long		$Torsional_J_inside_long		$Iy_inside_long		$Iz_inside_long_mod		$lineartag_Vertical
            var ele302 = new OpenSees.Elements.ElasticBeam3dWrapper(302, Area_inside_long, E_Timber_Long, G_Timber_Long, Torsional_J_inside_long, Iy_inside_long, Iz_inside_long_mod, 30801, 40501, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele302);

            //element elasticBeamColumn 303	30101		40101		 $Area_inside_long		$E_Timber_Long		$G_Timber_Long		$Torsional_J_inside_long		$Iy_inside_long		$Iz_inside_long_mod		$lineartag_Vertical
            var ele303 = new OpenSees.Elements.ElasticBeam3dWrapper(303, Area_inside_long, E_Timber_Long, G_Timber_Long, Torsional_J_inside_long, Iy_inside_long, Iz_inside_long_mod, 30101, 40101, lineartag_VerticalCrdT, 0, 0, 0);
            theDomain.AddElement(ele303);

            //element elasticBeamColumn	75	201		106		$Area_inside_trans_Rec	$E_Timber_Trans 	$G_Timber_Trans 	$Torsional_J_inside_trans_Rec		$Iy_inside_trans_Rec		$Iz_inside_trans_Rec_mod		$lineartag_Horizontal
            var ele75 = new OpenSees.Elements.ElasticBeam3dWrapper(75, Area_inside_trans_Rec, E_Timber_Trans, G_Timber_Trans, Torsional_J_inside_trans_Rec, Iy_inside_trans_Rec, Iz_inside_trans_Rec_mod, 201, 106, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele75);

            //element elasticBeamColumn	76	207		107		$Area_inside_trans_Rec	$E_Timber_Trans 	$G_Timber_Trans 	$Torsional_J_inside_trans_Rec		$Iy_inside_trans_Rec		$Iz_inside_trans_Rec_mod		$lineartag_Horizontal
            var ele76 = new OpenSees.Elements.ElasticBeam3dWrapper(76, Area_inside_trans_Rec, E_Timber_Trans, G_Timber_Trans, Torsional_J_inside_trans_Rec, Iy_inside_trans_Rec, Iz_inside_trans_Rec_mod, 207, 107, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele76);

            //element elasticBeamColumn	77	213		113		$Area_inside_trans_Rec	$E_Timber_Trans 	$G_Timber_Trans 	$Torsional_J_inside_trans_Rec		$Iy_inside_trans_Rec		$Iz_inside_trans_Rec_mod		$lineartag_Horizontal
            var ele77 = new OpenSees.Elements.ElasticBeam3dWrapper(77, Area_inside_trans_Rec, E_Timber_Trans, G_Timber_Trans, Torsional_J_inside_trans_Rec, Iy_inside_trans_Rec, Iz_inside_trans_Rec_mod, 213, 113, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele77);

            //element elasticBeamColumn	78	214		119		$Area_inside_trans_Rec	$E_Timber_Trans 	$G_Timber_Trans 	$Torsional_J_inside_trans_Rec		$Iy_inside_trans_Rec		$Iz_inside_trans_Rec_mod		$lineartag_Horizontal
            var ele78 = new OpenSees.Elements.ElasticBeam3dWrapper(78, Area_inside_trans_Rec, E_Timber_Trans, G_Timber_Trans, Torsional_J_inside_trans_Rec, Iy_inside_trans_Rec, Iz_inside_trans_Rec_mod, 214, 119, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele78);

            //element elasticBeamColumn	304	20101	10601	$Area_inside_trans_Rec	$E_Timber_Trans 	$G_Timber_Trans 	$Torsional_J_inside_trans_Rec		$Iy_inside_trans_Rec		$Iz_inside_trans_Rec_mod		$lineartag_Horizontal
            var ele304 = new OpenSees.Elements.ElasticBeam3dWrapper(304, Area_inside_trans_Rec, E_Timber_Trans, G_Timber_Trans, Torsional_J_inside_trans_Rec, Iy_inside_trans_Rec, Iz_inside_trans_Rec_mod, 20101, 10601, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele304);

            //element elasticBeamColumn	305	20701	10701	$Area_inside_trans_Rec	$E_Timber_Trans 	$G_Timber_Trans 	$Torsional_J_inside_trans_Rec		$Iy_inside_trans_Rec		$Iz_inside_trans_Rec_mod		$lineartag_Horizontal
            var ele305 = new OpenSees.Elements.ElasticBeam3dWrapper(305, Area_inside_trans_Rec, E_Timber_Trans, G_Timber_Trans, Torsional_J_inside_trans_Rec, Iy_inside_trans_Rec, Iz_inside_trans_Rec_mod, 20701, 10701, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele305);

            //element elasticBeamColumn	306	21301	11301	$Area_inside_trans_Rec	$E_Timber_Trans 	$G_Timber_Trans 	$Torsional_J_inside_trans_Rec		$Iy_inside_trans_Rec		$Iz_inside_trans_Rec_mod		$lineartag_Horizontal
            var ele306 = new OpenSees.Elements.ElasticBeam3dWrapper(306, Area_inside_trans_Rec, E_Timber_Trans, G_Timber_Trans, Torsional_J_inside_trans_Rec, Iy_inside_trans_Rec, Iz_inside_trans_Rec_mod, 21301, 11301, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele306);

            //element elasticBeamColumn	307	21401	11901	$Area_inside_trans_Rec	$E_Timber_Trans 	$G_Timber_Trans 	$Torsional_J_inside_trans_Rec		$Iy_inside_trans_Rec		$Iz_inside_trans_Rec_mod		$lineartag_Horizontal
            var ele307 = new OpenSees.Elements.ElasticBeam3dWrapper(307, Area_inside_trans_Rec, E_Timber_Trans, G_Timber_Trans, Torsional_J_inside_trans_Rec, Iy_inside_trans_Rec, Iz_inside_trans_Rec_mod, 21401, 11901, lineartag_HorizontalCrdT, 0, 0, 0);
            theDomain.AddElement(ele307);

            {
                //element twoNodeLink		79		2020		202		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele79 = new OpenSees.Elements.TwoNodeLinkWrapper(79, 3, 2020, 202, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele79);
            }
            {
                //element twoNodeLink		80		2030		203		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele80 = new OpenSees.Elements.TwoNodeLinkWrapper(80, 3, 2030, 203, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele80);
            }
            {
                //element twoNodeLink		81		2040		204		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele81 = new OpenSees.Elements.TwoNodeLinkWrapper(81, 3, 2040, 204, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele81);
            }
            {
                //element twoNodeLink		82		2050		205		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele82 = new OpenSees.Elements.TwoNodeLinkWrapper(82, 3, 2050, 205, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele82);
            }
            {
                //element twoNodeLink		83		2060		206		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele83 = new OpenSees.Elements.TwoNodeLinkWrapper(83, 3, 2060, 206, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele83);
            }
            {
                //element twoNodeLink		84		2080		208		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele84 = new OpenSees.Elements.TwoNodeLinkWrapper(84, 3, 2080, 208, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele84);
            }
            {
                //element twoNodeLink		85		2090		209		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele85 = new OpenSees.Elements.TwoNodeLinkWrapper(85, 3, 2090, 209, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele85);
            }
            {
                //element twoNodeLink		86		2100		210		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele86 = new OpenSees.Elements.TwoNodeLinkWrapper(86, 3, 2100, 210, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele86);
            }
            {
                //element twoNodeLink		87		2110		211		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele87 = new OpenSees.Elements.TwoNodeLinkWrapper(87, 3, 2110, 211, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele87);
            }
            {
                //element twoNodeLink		88		2120		212		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele88 = new OpenSees.Elements.TwoNodeLinkWrapper(88, 3, 2120, 212, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele88);
            }
            {
                //element twoNodeLink		89		2150		215		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele89 = new OpenSees.Elements.TwoNodeLinkWrapper(89, 3, 2150, 215, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele89);
            }
            {
                //element twoNodeLink		90		2160		216		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele90 = new OpenSees.Elements.TwoNodeLinkWrapper(90, 3, 2160, 216, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele90);
            }
            {
                //element twoNodeLink		91		2170		217		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele91 = new OpenSees.Elements.TwoNodeLinkWrapper(91, 3, 2170, 217, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele91);
            }
            {
                //element twoNodeLink		92		2180		218		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele92 = new OpenSees.Elements.TwoNodeLinkWrapper(92, 3, 2180, 218, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele92);
            }
            {
                //element twoNodeLink		93		2190		219		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele93 = new OpenSees.Elements.TwoNodeLinkWrapper(93, 3, 2190, 219, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele93);
            }
            {
                //element twoNodeLink		94		101			1010	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele94 = new OpenSees.Elements.TwoNodeLinkWrapper(94, 3, 101, 1010, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele94);
            }
            {
                //element twoNodeLink		95		102			1020	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele95 = new OpenSees.Elements.TwoNodeLinkWrapper(95, 3, 102, 1020, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele95);
            }
            {
                //element twoNodeLink		96		103			1030	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele96 = new OpenSees.Elements.TwoNodeLinkWrapper(96, 3, 103, 1030, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele96);
            }
            {
                //element twoNodeLink		97		104			1040	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele97 = new OpenSees.Elements.TwoNodeLinkWrapper(97, 3, 104, 1040, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele97);
            }
            {
                //element twoNodeLink		98		105			1050	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele98 = new OpenSees.Elements.TwoNodeLinkWrapper(98, 3, 105, 1050, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele98);
            }
            {
                //element twoNodeLink		99		108			1080	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele99 = new OpenSees.Elements.TwoNodeLinkWrapper(99, 3, 108, 1080, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele99);
            }
            {
                //element twoNodeLink		100		109			1090	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele100 = new OpenSees.Elements.TwoNodeLinkWrapper(100, 3, 109, 1090, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele100);
            }
            {
                //element twoNodeLink		101		110			1100	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele101 = new OpenSees.Elements.TwoNodeLinkWrapper(101, 3, 110, 1100, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele101);
            }
            {
                //element twoNodeLink		102		111			1110	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele102 = new OpenSees.Elements.TwoNodeLinkWrapper(102, 3, 111, 1110, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele102);
            }
            {
                //element twoNodeLink		103		112			1120	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele103 = new OpenSees.Elements.TwoNodeLinkWrapper(103, 3, 112, 1120, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele103);
            }
            {
                //element twoNodeLink		104		114			1140	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele104 = new OpenSees.Elements.TwoNodeLinkWrapper(104, 3, 114, 1140, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele104);
            }
            {
                //element twoNodeLink		105		115			1150	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele105 = new OpenSees.Elements.TwoNodeLinkWrapper(105, 3, 115, 1150, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele105);
            }
            {
                //element twoNodeLink		106		116			1160	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele106 = new OpenSees.Elements.TwoNodeLinkWrapper(106, 3, 116, 1160, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele106);
            }
            {
                //element twoNodeLink		107		117			1170	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele107 = new OpenSees.Elements.TwoNodeLinkWrapper(107, 3, 117, 1170, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele107);
            }
            {
                //element twoNodeLink		108		118			1180	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	0	0	1
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele108 = new OpenSees.Elements.TwoNodeLinkWrapper(108, 3, 118, 1180, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 0, 0, 1 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele108);
            }
            {
                //element twoNodeLink		109		3020		302		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	1	0	0
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele109 = new OpenSees.Elements.TwoNodeLinkWrapper(109, 3, 3020, 302, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 1, 0, 0 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele109);
            }
            {
                //element twoNodeLink		110		3040		304		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	1	0	0
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele110 = new OpenSees.Elements.TwoNodeLinkWrapper(110, 3, 3040, 304, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 1, 0, 0 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele110);
            }
            {
                //element twoNodeLink		111		3050		305		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	1	0	0
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele111 = new OpenSees.Elements.TwoNodeLinkWrapper(111, 3, 3050, 305, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 1, 0, 0 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele111);
            }
            {
                //element twoNodeLink		112		3060		306		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	1	0	0
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele112 = new OpenSees.Elements.TwoNodeLinkWrapper(112, 3, 3060, 306, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 1, 0, 0 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele112);
            }
            {
                //element twoNodeLink		113		3070		307		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	1	0	0
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele113 = new OpenSees.Elements.TwoNodeLinkWrapper(113, 3, 3070, 307, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 1, 0, 0 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele113);
            }
            {
                //element twoNodeLink		114		3090		309		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	1	0	0
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele114 = new OpenSees.Elements.TwoNodeLinkWrapper(114, 3, 3090, 309, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 1, 0, 0 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele114);
            }
            {
                //element twoNodeLink		115		3110		311		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	1	0	0
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele115 = new OpenSees.Elements.TwoNodeLinkWrapper(115, 3, 3110, 311, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 1, 0, 0 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele115);
            }
            {
                //element twoNodeLink		116		3120		312		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	1	0	0
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele116 = new OpenSees.Elements.TwoNodeLinkWrapper(116, 3, 3120, 312, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 1, 0, 0 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele116);
            }
            {
                //element twoNodeLink		117		3130		313		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	1	0	0
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele117 = new OpenSees.Elements.TwoNodeLinkWrapper(117, 3, 3130, 313, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 1, 0, 0 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele117);
            }
            {
                //element twoNodeLink		118		3140		314		-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	1	0	0
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele118 = new OpenSees.Elements.TwoNodeLinkWrapper(118, 3, 3140, 314, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 1, 0, 0 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele118);
            }
            {
                //element twoNodeLink		119		402			4020	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	1	0	0
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele119 = new OpenSees.Elements.TwoNodeLinkWrapper(119, 3, 402, 4020, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 1, 0, 0 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele119);
            }
            {
                //element twoNodeLink		120		403			4030	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	1	0	0
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele120 = new OpenSees.Elements.TwoNodeLinkWrapper(120, 3, 403, 4030, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 1, 0, 0 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele120);
            }
            {
                //element twoNodeLink		121		404			4040	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	1	0	0
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele121 = new OpenSees.Elements.TwoNodeLinkWrapper(121, 3, 404, 4040, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 1, 0, 0 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele121);
            }
            {
                //element twoNodeLink		122		406			4060	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	1	0	0
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele122 = new OpenSees.Elements.TwoNodeLinkWrapper(122, 3, 406, 4060, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 1, 0, 0 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele122);
            }
            {
                //element twoNodeLink		123		407			4070	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	1	0	0
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele123 = new OpenSees.Elements.TwoNodeLinkWrapper(123, 3, 407, 4070, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 1, 0, 0 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele123);
            }
            {
                //element twoNodeLink		124		409			4090	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	1	0	0
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele124 = new OpenSees.Elements.TwoNodeLinkWrapper(124, 3, 409, 4090, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 1, 0, 0 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele124);
            }
            {
                //element twoNodeLink		125		411			4110	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	1	0	0
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele125 = new OpenSees.Elements.TwoNodeLinkWrapper(125, 3, 411, 4110, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 1, 0, 0 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele125);
            }
            {
                //element twoNodeLink		126		413			4130	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	1	0	0
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele126 = new OpenSees.Elements.TwoNodeLinkWrapper(126, 3, 413, 4130, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 1, 0, 0 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele126);
            }
            {
                //element twoNodeLink		127		414			4140	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	1	0	0
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele127 = new OpenSees.Elements.TwoNodeLinkWrapper(127, 3, 414, 4140, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 1, 0, 0 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele127);
            }
            {
                //element twoNodeLink		128		412			4120	-mat	$Stiffness_Tension	$Stiffness_IP	$Stiffness_OOP	$Stiffness_about_Tension	$Stiffness_about_OOP	$Stiffness_about_IP	-dir	1	2	3	4	5	6 -orient	1	0	0
                var jointElementMatsSet = new OpenSees.Materials.Uniaxials.UniaxialMaterialWrapper[] { Stiffness_TensionMat, Stiffness_IPMat, Stiffness_OOPMat, Stiffness_about_TensionMat, Stiffness_about_OOPMat, Stiffness_about_IPMat };
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, 1, 2, 3, 4, 5 });
                var tlele128 = new OpenSees.Elements.TwoNodeLinkWrapper(128, 3, 412, 4120, dirIds, jointElementMatsSet, new OpenSees.VectorWrapper(new double[] { 1, 0, 0 }), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), new OpenSees.VectorWrapper(0), 0, 0);
                theDomain.AddElement(tlele128);
            }

            #endregion

            #region recorders
            var savepath = System.Environment.CurrentDirectory + @"\opsnet_results\Displacements\";
            if (!System.IO.Directory.Exists(savepath))
                System.IO.Directory.CreateDirectory(savepath);

            {
                var dirIds = new OpenSees.IDWrapper(new int[] { 0, });
                var nodeTags = new int[] { 412,404,217,210,204 };
                var counter = 1;
                foreach (var tag in nodeTags)
                {
                    var opsstream = new OpenSees.Handlers.DataFileStreamWrapper(savepath + $@"\RF{counter++}.out");
                    var recorder = new OpenSees.Recorders.NodeRecorderWrapper(dirIds, new OpenSees.IDWrapper(new int[] { tag }), 0, "disp", theDomain, opsstream);
                    theDomain.AddRecorder(recorder);
                }
            }

            #endregion

            #region loading 
            var theSeries = new OpenSees.Components.Timeseries.LinearSeriesWrapper();
            var theLoadPattern = new OpenSees.Components.LoadPatterns.LoadPatternWrapper(1);
            theLoadPattern.SetTimeSeries(theSeries);
            theDomain.AddLoadPattern(theLoadPattern);

            var theLoadValues = new OpenSees.VectorWrapper(6);
            theLoadValues.Zero();
            theLoadValues[0] = 920.992;

            theDomain.AddNodalLoad(new OpenSees.Components.Loads.NodalLoadWrapper(1, 2020, theLoadValues, false), 1);
            theDomain.AddNodalLoad(new OpenSees.Components.Loads.NodalLoadWrapper(2, 2030, theLoadValues, false), 1);
            theDomain.AddNodalLoad(new OpenSees.Components.Loads.NodalLoadWrapper(3, 2040, theLoadValues, false), 1);
            theDomain.AddNodalLoad(new OpenSees.Components.Loads.NodalLoadWrapper(4, 2050, theLoadValues, false), 1);
            theDomain.AddNodalLoad(new OpenSees.Components.Loads.NodalLoadWrapper(5, 2060, theLoadValues, false), 1);
            theDomain.AddNodalLoad(new OpenSees.Components.Loads.NodalLoadWrapper(6, 2080, theLoadValues, false), 1);
            theDomain.AddNodalLoad(new OpenSees.Components.Loads.NodalLoadWrapper(7, 2090, theLoadValues, false), 1);
            theDomain.AddNodalLoad(new OpenSees.Components.Loads.NodalLoadWrapper(8, 2100, theLoadValues, false), 1);
            theDomain.AddNodalLoad(new OpenSees.Components.Loads.NodalLoadWrapper(9, 2110, theLoadValues, false), 1);
            theDomain.AddNodalLoad(new OpenSees.Components.Loads.NodalLoadWrapper(10, 2120, theLoadValues, false), 1);
            theDomain.AddNodalLoad(new OpenSees.Components.Loads.NodalLoadWrapper(11, 2150, theLoadValues, false), 1);
            theDomain.AddNodalLoad(new OpenSees.Components.Loads.NodalLoadWrapper(12, 2160, theLoadValues, false), 1);
            theDomain.AddNodalLoad(new OpenSees.Components.Loads.NodalLoadWrapper(13, 2170, theLoadValues, false), 1);
            theDomain.AddNodalLoad(new OpenSees.Components.Loads.NodalLoadWrapper(14, 2180, theLoadValues, false), 1);
            theDomain.AddNodalLoad(new OpenSees.Components.Loads.NodalLoadWrapper(15, 2190, theLoadValues, false), 1);
            theDomain.AddNodalLoad(new OpenSees.Components.Loads.NodalLoadWrapper(16, 4020, theLoadValues, false), 1);
            theDomain.AddNodalLoad(new OpenSees.Components.Loads.NodalLoadWrapper(17, 4030, theLoadValues, false), 1);
            theDomain.AddNodalLoad(new OpenSees.Components.Loads.NodalLoadWrapper(18, 4040, theLoadValues, false), 1);
            theDomain.AddNodalLoad(new OpenSees.Components.Loads.NodalLoadWrapper(19, 4060, theLoadValues, false), 1);
            theDomain.AddNodalLoad(new OpenSees.Components.Loads.NodalLoadWrapper(20, 4070, theLoadValues, false), 1);
            theDomain.AddNodalLoad(new OpenSees.Components.Loads.NodalLoadWrapper(21, 4090, theLoadValues, false), 1);
            theDomain.AddNodalLoad(new OpenSees.Components.Loads.NodalLoadWrapper(22, 4110, theLoadValues, false), 1);
            theDomain.AddNodalLoad(new OpenSees.Components.Loads.NodalLoadWrapper(23, 4130, theLoadValues, false), 1);
            theDomain.AddNodalLoad(new OpenSees.Components.Loads.NodalLoadWrapper(24, 4140, theLoadValues, false), 1);
            theDomain.AddNodalLoad(new OpenSees.Components.Loads.NodalLoadWrapper(25, 4120, theLoadValues, false), 1);


            var theModel = new OpenSees.AnalysisModelWrapper();
            var theSolnAlgo = new OpenSees.Algorithms.NewtonRaphsonWrapper();
            var theIntegrator = new OpenSees.Integrators.Static.LoadControlWrapper(0.05, 1, 0.05, 0.05);
            var theHandler = new OpenSees.Handlers.TransformationConstraintHandlerWrapper();
            var theRCM = new OpenSees.GraphNumberers.RCMWrapper(false);
            var theNumberer = new OpenSees.Numberers.DOF_NumbererWrapper(theRCM);
            var theSolver = new OpenSees.Systems.Linears.BandGenLinLapackSolverWrapper();
            var theSOE = new OpenSees.Systems.Linears.BandGenLinSOEWrapper(theSolver);
            var test = new OpenSees.ConvergenceTests.CTestNormDispIncrWrapper(1e-6, 1000, 2, 2, 1.0e10);
            var theAnalysis = new OpenSees.Analysis.StaticAnalysisWrapper(
                theDomain,
                theHandler,
                theNumberer,
                theModel,
                theSolnAlgo,
                theSOE,
                theIntegrator,
                test);

            theAnalysis.Analyze(20);


            Console.ReadKey();
            #endregion
        }