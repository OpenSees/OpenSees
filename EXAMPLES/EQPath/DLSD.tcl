puts " REF:"
puts " Structural Engineering and Mechanics, Vol. 48, No. 6 (2013) 879-914"
puts " DOI: http://dx.doi.org/10.12989/sem.2013.48.6.879"
puts " Comprehensive evaluation of structural geometrical nonlinear solution techniques "
puts " Part II: Comparing efficiencies of the methods"
puts " M. Rezaiee-Pajand , M. Ghalishooyan and M. Salehi-Ahmadabad "
puts " 3.14 Problem fourteen page 27 "
puts " double layer circular dome"

wipe

model basic -ndm 3 -ndf 3;

node	0	13	0	0
node	1	12.71591889	2.702851596	0
node	2	11.87609127	5.287575642	0
node	3	10.51722162	7.641207326	0
node	4	8.698699051	9.660881679	0
node	5	6.500001702	11.25832927	0
node	6	4.01722317	12.36373398	0
node	7	1.358872759	12.92878435	0
node	8	-1.358866895	12.92878497	0
node	9	-4.017217562	12.36373581	0
node	10	-6.499996596	11.25833221	0
node	11	-8.698694669	9.660885625	0
node	12	-10.51721815	7.641212096	0
node	13	-11.87608887	5.287581029	0
node	14	-12.71591767	2.702857364	0
node	15	-13	0	0
node	16	-12.71592012	-2.702845828	0
node	17	-11.87609367	-5.287570255	0
node	18	-10.51722509	-7.641202555	0
node	19	-8.698703433	-9.660877733	0
node	20	-6.500006809	-11.25832632	0
node	21	-4.017228778	-12.36373216	0
node	22	-1.358878624	-12.92878374	0
node	23	1.35886103	-12.92878558	0
node	24	4.017211954	-12.36373763	0
node	25	6.499991489	-11.25833516	0
node	26	8.698690287	-9.66088957	0
node	27	10.51721469	-7.641216867	0
node	28	11.87608647	-5.287586415	0
node	29	12.71591644	-2.702863132	0
node	30	11.93426276	1.254341379	-0.51
node	31	11.41267836	3.708203415	-0.51
node	32	10.3923053	5.999999214	-0.51
node	33	8.917738756	8.029566332	-0.51
node	34	7.053424349	9.708202973	-0.51
node	35	4.88084154	10.96254468	-0.51
node	36	2.494942597	11.73777072	-0.51
node	37	0	12	-0.51
node	38	-2.494937273	11.73777185	-0.51
node	39	-4.880836568	10.96254689	-0.51
node	40	-7.053419945	9.708206172	-0.51
node	41	-8.917735113	8.029570377	-0.51
node	42	-10.39230258	6.000003928	-0.51
node	43	-11.41267668	3.708208592	-0.51
node	44	-11.93426219	1.254346792	-0.51
node	45	-11.93426333	-1.254335966	-0.51
node	46	-11.41268005	-3.708198238	-0.51
node	47	-10.39230802	-5.999994501	-0.51
node	48	-8.917742398	-8.029562287	-0.51
node	49	-7.053428752	-9.708199773	-0.51
node	50	-4.880846513	-10.96254247	-0.51
node	51	-2.494947921	-11.73776959	-0.51
node	52	0	-12	-0.51
node	53	2.494931949	-11.73777298	-0.51
node	54	4.880831595	-10.96254911	-0.51
node	55	7.053415541	-9.708209371	-0.51
node	56	8.917731471	-8.029574422	-0.51
node	57	10.39229986	-6.000008642	-0.51
node	58	11.412675	-3.708213768	-0.51
node	59	11.93426163	-1.254352205	-0.51
node	60	11	0	1
node	61	10.75962368	2.287028274	1
node	62	10.0490003	4.474102466	1
node	63	8.899187525	6.465636968	1
node	64	7.360437659	8.17459219	1
node	65	5.50000144	9.52627861	1
node	66	3.399188836	10.46162106	1
node	67	1.149815412	10.93974061	1
node	68	-1.149810449	10.93974113	1
node	69	-3.399184091	10.4616226	1
node	70	-5.499997119	9.526281105	1
node	71	-7.360433951	8.174595529	1
node	72	-8.899184592	6.465641004	1
node	73	-10.04899828	4.474107024	1
node	74	-10.75962264	2.287033154	1
node	75	-11	0	1
node	76	-10.75962471	-2.287023393	1
node	77	-10.04900233	-4.474097908	1
node	78	-8.899190457	-6.465632931	1
node	79	-7.360441367	-8.174588851	1
node	80	-5.500005761	-9.526276115	1
node	81	-3.399193582	-10.46161952	1
node	82	-1.149820374	-10.93974008	1
node	83	1.149805487	-10.93974165	1
node	84	3.399179346	-10.46162415	1
node	85	5.499992798	-9.5262836	1
node	86	7.360430243	-8.174598867	1
node	87	8.899181659	-6.465645041	1
node	88	10.04899625	-4.474111582	1
node	89	10.7596216	-2.287038035	1
node	90	9.945218969	1.045284482	0.37
node	91	9.510565303	3.090169512	0.37
node	92	8.660254416	4.999999345	0.37
node	93	7.431448963	6.691305277	0.37
node	94	5.877853624	8.090169144	0.37
node	95	4.06736795	9.1354539	0.37
node	96	2.079118831	9.781475599	0.37
node	97	0	10	0.37
node	98	-2.079114394	9.781476542	0.37
node	99	-4.067363806	9.135455745	0.37
node	100	-5.877849954	8.09017181	0.37
node	101	-7.431445928	6.691308648	0.37
node	102	-8.660252148	5.000003274	0.37
node	103	-9.510563901	3.090173826	0.37
node	104	-9.945218495	1.045288993	0.37
node	105	-9.945219444	-1.045279971	0.37
node	106	-9.510566705	-3.090165198	0.37
node	107	-8.660256684	-4.999995417	0.37
node	108	-7.431451998	-6.691301906	0.37
node	109	-5.877857293	-8.090166478	0.37
node	110	-4.067372094	-9.135452055	0.37
node	111	-2.079123268	-9.781474656	0.37
node	112	0	-10	0.37
node	113	2.079109957	-9.781477485	0.37
node	114	4.067359663	-9.13545759	0.37
node	115	5.877846285	-8.090174476	0.37
node	116	7.431442893	-6.691312019	0.37
node	117	8.66024988	-5.000007202	0.37
node	118	9.5105625	-3.09017814	0.37
node	119	9.945218021	-1.045293504	0.37
node	120	9	0	1.75
node	121	8.803328463	1.871204951	1.75
node	122	8.22190934	3.66062929	1.75
node	123	7.281153429	5.29006661	1.75
node	124	6.022176266	6.688302701	1.75
node	125	4.500001178	7.794227954	1.75
node	126	2.781154502	8.559508142	1.75
node	127	0.940758064	8.950696859	1.75
node	128	-0.940754004	8.950697286	1.75
node	129	-2.78115062	8.559509404	1.75
node	130	-4.499997643	7.794229995	1.75
node	131	-6.022173232	6.688305432	1.75
node	132	-7.28115103	5.290069913	1.75
node	133	-8.22190768	3.66063302	1.75
node	134	-8.803327614	1.871208944	1.75
node	135	-9	0	1.75
node	136	-8.803329312	-1.871200958	1.75
node	137	-8.221911001	-3.660625561	1.75
node	138	-7.281155829	-5.290063307	1.75
node	139	-6.0221793	-6.688299969	1.75
node	140	-4.500004714	-7.794225913	1.75
node	141	-2.781158385	-8.559506881	1.75
node	142	-0.940762124	-8.950696432	1.75
node	143	0.940749944	-8.950697713	1.75
node	144	2.781146737	-8.559510665	1.75
node	145	4.499994108	-7.794232036	1.75
node	146	6.022170199	-6.688308164	1.75
node	147	7.28114863	-5.290073215	1.75
node	148	8.221906019	-3.660636749	1.75
node	149	8.803326766	-1.871212937	1.75
node	150	7.956175176	0.836227586	1.12
node	151	7.608452242	2.47213561	1.12
node	152	6.928203533	3.999999476	1.12
node	153	5.94515917	5.353044222	1.12
node	154	4.702282899	6.472135315	1.12
node	155	3.25389436	7.30836312	1.12
node	156	1.663295065	7.825180479	1.12
node	157	0	8	1.12
node	158	-1.663291515	7.825181233	1.12
node	159	-3.253891045	7.308364596	1.12
node	160	-4.702279963	6.472137448	1.12
node	161	-5.945156742	5.353046918	1.12
node	162	-6.928201718	4.000002619	1.12
node	163	-7.608451121	2.472139061	1.12
node	164	-7.956174796	0.836231195	1.12
node	165	-7.956175555	-0.836223977	1.12
node	166	-7.608453364	-2.472132159	1.12
node	167	-6.928205347	-3.999996334	1.12
node	168	-5.945161598	-5.353041525	1.12
node	169	-4.702285835	-6.472133182	1.12
node	170	-3.253897675	-7.308361644	1.12
node	171	-1.663298614	-7.825179724	1.12
node	172	0	-8	1.12
node	173	1.663287966	-7.825181988	1.12
node	174	3.25388773	-7.308366072	1.12
node	175	4.702277028	-6.472139581	1.12
node	176	5.945154314	-5.353049615	1.12
node	177	6.928199904	-4.000005761	1.12
node	178	7.60845	-2.472142512	1.12
node	179	7.956174417	-0.836234804	1.12
node	180	7	0	2.5
node	181	6.847033249	1.455381629	2.5
node	182	6.394818376	2.847156115	2.5
node	183	5.663119334	4.114496252	2.5
node	184	4.683914874	5.202013212	2.5
node	185	3.500000917	6.062177297	2.5
node	186	2.163120169	6.657395222	2.5
node	187	0.731700716	6.961653113	2.5
node	188	-0.731697559	6.961653445	2.5
node	189	-2.163117149	6.657396203	2.5
node	190	-3.499998167	6.062178885	2.5
node	191	-4.683912514	5.202015336	2.5
node	192	-5.663117468	4.114498821	2.5
node	193	-6.394817084	2.847159015	2.5
node	194	-6.847032589	1.455384734	2.5
node	195	-7	0	2.5
node	196	-6.847033909	-1.455378523	2.5
node	197	-6.394819667	-2.847153214	2.5
node	198	-5.6631212	-4.114493684	2.5
node	199	-4.683917233	-5.202011087	2.5
node	200	-3.500003666	-6.06217571	2.5
node	201	-2.163123188	-6.65739424	2.5
node	202	-0.731703874	-6.961652781	2.5
node	203	0.731694401	-6.961653776	2.5
node	204	2.163114129	-6.657397184	2.5
node	205	3.499995417	-6.062180472	2.5
node	206	4.683910155	-5.202017461	2.5
node	207	5.663115601	-4.11450139	2.5
node	208	6.394815793	-2.847161916	2.5
node	209	6.847031929	-1.45538784	2.5
node	210	5.967131382	0.627170689	1.74
node	211	5.706339182	1.854101707	1.74
node	212	5.19615265	2.999999607	1.74
node	213	4.458869378	4.014783166	1.74
node	214	3.526712174	4.854101486	1.74
node	215	2.44042077	5.48127234	1.74
node	216	1.247471298	5.868885359	1.74
node	217	0	6	1.74
node	218	-1.247468636	5.868885925	1.74
node	219	-2.440418284	5.481273447	1.74
node	220	-3.526709973	4.854103086	1.74
node	221	-4.458867557	4.014785189	1.74
node	222	-5.196151289	3.000001964	1.74
node	223	-5.706338341	1.854104296	1.74
node	224	-5.967131097	0.627173396	1.74
node	225	-5.967131666	-0.627167983	1.74
node	226	-5.706340023	-1.854099119	1.74
node	227	-5.19615401	-2.99999725	1.74
node	228	-4.458871199	-4.014781144	1.74
node	229	-3.526714376	-4.854099887	1.74
node	230	-2.440423256	-5.481271233	1.74
node	231	-1.247473961	-5.868884793	1.74
node	232	0	-6	1.74
node	233	1.247465974	-5.868886491	1.74
node	234	2.440415798	-5.481274554	1.74
node	235	3.526707771	-4.854104686	1.74
node	236	4.458865736	-4.014787211	1.74
node	237	5.196149928	-3.000004321	1.74
node	238	5.7063375	-1.854106884	1.74
node	239	5.967130813	-0.627176103	1.74
node	240	5	0	3
node	241	4.890738035	1.039558306	3
node	242	4.567727411	2.033682939	3
node	243	4.045085238	2.938925895	3
node	244	3.345653481	3.715723723	3
node	245	2.500000655	4.330126641	3
node	246	1.545085835	4.755282301	3
node	247	0.522643369	4.972609366	3
node	248	-0.522641113	4.972609603	3
node	249	-1.545083678	4.755283002	3
node	250	-2.499998691	4.330127775	3
node	251	-3.345651796	3.71572524	3
node	252	-4.045083905	2.938927729	3
node	253	-4.567726489	2.033685011	3
node	254	-4.890737564	1.039560525	3
node	255	-5	0	3
node	256	-4.890738507	-1.039556088	3
node	257	-4.567728334	-2.033680867	3
node	258	-4.045086572	-2.93892406	3
node	259	-3.345655167	-3.715722205	3
node	260	-2.500002619	-4.330125507	3
node	261	-1.545087992	-4.7552816	3
node	262	-0.522645624	-4.972609129	3
node	263	0.522638858	-4.97260984	3
node	264	1.545081521	-4.755283703	3
node	265	2.499996726	-4.330128909	3
node	266	3.34565011	-3.715726758	3
node	267	4.045082572	-2.938929564	3
node	268	4.567725566	-2.033687083	3
node	269	4.890737092	-1.039562743	3
node	270	3.978087588	0.418113793	2.42
node	271	3.804226121	1.236067805	2.42
node	272	3.464101766	1.999999738	2.42
node	273	2.972579585	2.676522111	2.42
node	274	2.35114145	3.236067658	2.42
node	275	1.62694718	3.65418156	2.42
node	276	0.831647532	3.912590239	2.42
node	277	0	4	2.42
node	278	-0.831645758	3.912590617	2.42
node	279	-1.626945523	3.654182298	2.42
node	280	-2.351139982	3.236068724	2.42
node	281	-2.972578371	2.676523459	2.42
node	282	-3.464100859	2.000001309	2.42
node	283	-3.804225561	1.236069531	2.42
node	284	-3.978087398	0.418115597	2.42
node	285	-3.978087777	-0.418111989	2.42
node	286	-3.804226682	-1.236066079	2.42
node	287	-3.464102674	-1.999998167	2.42
node	288	-2.972580799	-2.676520762	2.42
node	289	-2.351142917	-3.236066591	2.42
node	290	-1.626948838	-3.654180822	2.42
node	291	-0.831649307	-3.912589862	2.42
node	292	0	-4	2.42
node	293	0.831643983	-3.912590994	2.42
node	294	1.626943865	-3.654183036	2.42
node	295	2.351138514	-3.23606979	2.42
node	296	2.972577157	-2.676524807	2.42
node	297	3.464099952	-2.000002881	2.42
node	298	3.804225	-1.236071256	2.42
node	299	3.978087208	-0.418117402	2.42
node	300	3	0	3.5
node	301	2.934442821	0.623734984	3.5
node	302	2.740636447	1.220209763	3.5
node	303	2.427051143	1.763355537	3.5
node	304	2.007392089	2.229434234	3.5
node	305	1.500000393	2.598075985	3.5
node	306	0.927051501	2.853169381	3.5
node	307	0.313586021	2.98356562	3.5
node	308	-0.313584668	2.983565762	3.5
node	309	-0.927050207	2.853169801	3.5
node	310	-1.499999214	2.598076665	3.5
node	311	-2.007391077	2.229435144	3.5
node	312	-2.427050343	1.763356638	3.5
node	313	-2.740635893	1.220211007	3.5
node	314	-2.934442538	0.623736315	3.5
node	315	-3	1.36E-06	3.5
node	316	-2.934443104	-0.623733653	3.5
node	317	-2.740637	-1.22020852	3.5
node	318	-2.427051943	-1.763354436	3.5
node	319	-2.0073931	-2.229433323	3.5
node	320	-1.500001571	-2.598075304	3.5
node	321	-0.927052795	-2.85316896	3.5
node	322	-0.313587375	-2.983565477	3.5
node	323	0.313583315	-2.983565904	3.5
node	324	0.927048912	-2.853170222	3.5
node	325	1.499998036	-2.598077345	3.5
node	326	2.007390066	-2.229436055	3.5
node	327	2.427049543	-1.763357738	3.5
node	328	2.74063534	-1.22021225	3.5
node	329	2.934442255	-0.623737646	3.5
node	330	1.989043794	0.209056896	2.62
node	331	1.902113061	0.618033902	2.62
node	332	1.732050883	0.999999869	2.62
node	333	1.486289793	1.338261055	2.62
node	334	1.175570725	1.618033829	2.62
node	335	0.81347359	1.82709078	2.62
node	336	0.415823766	1.95629512	2.62
node	337	0	2	2.62
node	338	-0.415822879	1.956295308	2.62
node	339	-0.813472761	1.827091149	2.62
node	340	-1.175569991	1.618034362	2.62
node	341	-1.486289186	1.33826173	2.62
node	342	-1.73205043	1.000000655	2.62
node	343	-1.90211278	0.618034765	2.62
node	344	-1.989043699	0.209057799	2.62
node	345	-1.989043889	-0.209055994	2.62
node	346	-1.902113341	-0.61803304	2.62
node	347	-1.732051337	-0.999999083	2.62
node	348	-1.4862904	-1.338260381	2.62
node	349	-1.175571459	-1.618033296	2.62
node	350	-0.813474419	-1.827090411	2.62
node	351	-0.415824654	-1.956294931	2.62
node	352	0	-2	2.62
node	353	0.415821991	-1.956295497	2.62
node	354	0.813471933	-1.827091518	2.62
node	355	1.175569257	-1.618034895	2.62
node	356	1.486288579	-1.338262404	2.62
node	357	1.732049976	-1.00000144	2.62
node	358	1.9021125	-0.618035628	2.62
node	359	1.989043604	-0.209058701	2.62
node	360	1	0	3.75
node	361	0.978147607	0.207911661	3.75
node	362	0.913545482	0.406736588	3.75
node	363	0.809017048	0.587785179	3.75
node	364	0.669130696	0.743144745	3.75
node	365	0.500000131	0.866025328	3.75
node	366	0.309017167	0.95105646	3.75
node	367	0.104528674	0.994521873	3.75
node	368	-0.104528223	0.994521921	3.75
node	369	-0.309016736	0.9510566	3.75
node	370	-0.499999738	0.866025555	3.75
node	371	-0.669130359	0.743145048	3.75
node	372	-0.809016781	0.587785546	3.75
node	373	-0.913545298	0.406737002	3.75
node	374	-0.978147513	0.207912105	3.75
node	375	-1	0	3.75
node	376	-0.978147701	-0.207911218	3.75
node	377	-0.913545667	-0.406736173	3.75
node	378	-0.809017314	-0.587784812	3.75
node	379	-0.669131033	-0.743144441	3.75
node	380	-0.500000524	-0.866025101	3.75
node	381	-0.309017598	-0.95105632	3.75
node	382	-0.104529125	-0.994521826	3.75
node	383	0.104527772	-0.994521968	3.75
node	384	0.309016304	-0.951056741	3.75
node	385	0.499999345	-0.866025782	3.75
node	386	0.669130022	-0.743145352	3.75
node	387	0.809016514	-0.587785913	3.75
node	388	0.913545113	-0.406737417	3.75
node	389	0.978147418	-0.207912549	3.75
node	390	0	0	3.88

set i 0
while {$i<30} {
    fix $i 1 1 1
	incr i
}
	

set A 450
set matTag 1
uniaxialMaterial Elastic $matTag 2.1e5

 element	corotTruss	0	30	0	$A	$matTag
element	corotTruss	1	31	1	$A	$matTag
element	corotTruss	2	32	2	$A	$matTag
element	corotTruss	3	33	3	$A	$matTag
element	corotTruss	4	34	4	$A	$matTag
element	corotTruss	5	35	5	$A	$matTag
element	corotTruss	6	36	6	$A	$matTag
element	corotTruss	7	37	7	$A	$matTag
element	corotTruss	8	38	8	$A	$matTag
element	corotTruss	9	39	9	$A	$matTag
element	corotTruss	10	40	10	$A	$matTag
element	corotTruss	11	41	11	$A	$matTag
element	corotTruss	12	42	12	$A	$matTag
element	corotTruss	13	43	13	$A	$matTag
element	corotTruss	14	44	14	$A	$matTag
element	corotTruss	15	45	15	$A	$matTag
element	corotTruss	16	46	16	$A	$matTag
element	corotTruss	17	47	17	$A	$matTag
element	corotTruss	18	48	18	$A	$matTag
element	corotTruss	19	49	19	$A	$matTag
element	corotTruss	20	50	20	$A	$matTag
element	corotTruss	21	51	21	$A	$matTag
element	corotTruss	22	52	22	$A	$matTag
element	corotTruss	23	53	23	$A	$matTag
element	corotTruss	24	54	24	$A	$matTag
element	corotTruss	25	55	25	$A	$matTag
element	corotTruss	26	56	26	$A	$matTag
element	corotTruss	27	57	27	$A	$matTag
element	corotTruss	28	58	28	$A	$matTag
element	corotTruss	29	59	29	$A	$matTag
element	corotTruss	30	1	30	$A	$matTag
element	corotTruss	31	2	31	$A	$matTag
element	corotTruss	32	3	32	$A	$matTag
element	corotTruss	33	4	33	$A	$matTag
element	corotTruss	34	5	34	$A	$matTag
element	corotTruss	35	6	35	$A	$matTag
element	corotTruss	36	7	36	$A	$matTag
element	corotTruss	37	8	37	$A	$matTag
element	corotTruss	38	9	38	$A	$matTag
element	corotTruss	39	10	39	$A	$matTag
element	corotTruss	40	11	40	$A	$matTag
element	corotTruss	41	12	41	$A	$matTag
element	corotTruss	42	13	42	$A	$matTag
element	corotTruss	43	14	43	$A	$matTag
element	corotTruss	44	15	44	$A	$matTag
element	corotTruss	45	16	45	$A	$matTag
element	corotTruss	46	17	46	$A	$matTag
element	corotTruss	47	18	47	$A	$matTag
element	corotTruss	48	19	48	$A	$matTag
element	corotTruss	49	20	49	$A	$matTag
element	corotTruss	50	21	50	$A	$matTag
element	corotTruss	51	22	51	$A	$matTag
element	corotTruss	52	23	52	$A	$matTag
element	corotTruss	53	24	53	$A	$matTag
element	corotTruss	54	25	54	$A	$matTag
element	corotTruss	55	26	55	$A	$matTag
element	corotTruss	56	27	56	$A	$matTag
element	corotTruss	57	28	57	$A	$matTag
element	corotTruss	58	29	58	$A	$matTag
element	corotTruss	59	0	59	$A	$matTag
element	corotTruss	60	31	30	$A	$matTag
element	corotTruss	61	32	31	$A	$matTag
element	corotTruss	62	33	32	$A	$matTag
element	corotTruss	63	34	33	$A	$matTag
element	corotTruss	64	35	34	$A	$matTag
element	corotTruss	65	36	35	$A	$matTag
element	corotTruss	66	37	36	$A	$matTag
element	corotTruss	67	38	37	$A	$matTag
element	corotTruss	68	39	38	$A	$matTag
element	corotTruss	69	40	39	$A	$matTag
element	corotTruss	70	41	40	$A	$matTag
element	corotTruss	71	42	41	$A	$matTag
element	corotTruss	72	43	42	$A	$matTag
element	corotTruss	73	44	43	$A	$matTag
element	corotTruss	74	45	44	$A	$matTag
element	corotTruss	75	46	45	$A	$matTag
element	corotTruss	76	47	46	$A	$matTag
element	corotTruss	77	48	47	$A	$matTag
element	corotTruss	78	49	48	$A	$matTag
element	corotTruss	79	50	49	$A	$matTag
element	corotTruss	80	51	50	$A	$matTag
element	corotTruss	81	52	51	$A	$matTag
element	corotTruss	82	53	52	$A	$matTag
element	corotTruss	83	54	53	$A	$matTag
element	corotTruss	84	55	54	$A	$matTag
element	corotTruss	85	56	55	$A	$matTag
element	corotTruss	86	57	56	$A	$matTag
element	corotTruss	87	58	57	$A	$matTag
element	corotTruss	88	59	58	$A	$matTag
element	corotTruss	89	30	59	$A	$matTag
element	corotTruss	90	60	0	$A	$matTag
element	corotTruss	91	61	1	$A	$matTag
element	corotTruss	92	62	2	$A	$matTag
element	corotTruss	93	63	3	$A	$matTag
element	corotTruss	94	64	4	$A	$matTag
element	corotTruss	95	65	5	$A	$matTag
element	corotTruss	96	66	6	$A	$matTag
element	corotTruss	97	67	7	$A	$matTag
element	corotTruss	98	68	8	$A	$matTag
element	corotTruss	99	69	9	$A	$matTag
element	corotTruss	100	70	10	$A	$matTag
element	corotTruss	101	71	11	$A	$matTag
element	corotTruss	102	72	12	$A	$matTag
element	corotTruss	103	73	13	$A	$matTag
element	corotTruss	104	74	14	$A	$matTag
element	corotTruss	105	75	15	$A	$matTag
element	corotTruss	106	76	16	$A	$matTag
element	corotTruss	107	77	17	$A	$matTag
element	corotTruss	108	78	18	$A	$matTag
element	corotTruss	109	79	19	$A	$matTag
element	corotTruss	110	80	20	$A	$matTag
element	corotTruss	111	81	21	$A	$matTag
element	corotTruss	112	82	22	$A	$matTag
element	corotTruss	113	83	23	$A	$matTag
element	corotTruss	114	84	24	$A	$matTag
element	corotTruss	115	85	25	$A	$matTag
element	corotTruss	116	86	26	$A	$matTag
element	corotTruss	117	87	27	$A	$matTag
element	corotTruss	118	88	28	$A	$matTag
element	corotTruss	119	89	29	$A	$matTag
element	corotTruss	120	60	30	$A	$matTag
element	corotTruss	121	61	31	$A	$matTag
element	corotTruss	122	62	32	$A	$matTag
element	corotTruss	123	63	33	$A	$matTag
element	corotTruss	124	64	34	$A	$matTag
element	corotTruss	125	65	35	$A	$matTag
element	corotTruss	126	66	36	$A	$matTag
element	corotTruss	127	67	37	$A	$matTag
element	corotTruss	128	68	38	$A	$matTag
element	corotTruss	129	69	39	$A	$matTag
element	corotTruss	130	70	40	$A	$matTag
element	corotTruss	131	71	41	$A	$matTag
element	corotTruss	132	72	42	$A	$matTag
element	corotTruss	133	73	43	$A	$matTag
element	corotTruss	134	74	44	$A	$matTag
element	corotTruss	135	75	45	$A	$matTag
element	corotTruss	136	76	46	$A	$matTag
element	corotTruss	137	77	47	$A	$matTag
element	corotTruss	138	78	48	$A	$matTag
element	corotTruss	139	79	49	$A	$matTag
element	corotTruss	140	80	50	$A	$matTag
element	corotTruss	141	81	51	$A	$matTag
element	corotTruss	142	82	52	$A	$matTag
element	corotTruss	143	83	53	$A	$matTag
element	corotTruss	144	84	54	$A	$matTag
element	corotTruss	145	85	55	$A	$matTag
element	corotTruss	146	86	56	$A	$matTag
element	corotTruss	147	87	57	$A	$matTag
element	corotTruss	148	88	58	$A	$matTag
element	corotTruss	149	89	59	$A	$matTag
element	corotTruss	150	90	30	$A	$matTag
element	corotTruss	151	91	31	$A	$matTag
element	corotTruss	152	92	32	$A	$matTag
element	corotTruss	153	93	33	$A	$matTag
element	corotTruss	154	94	34	$A	$matTag
element	corotTruss	155	95	35	$A	$matTag
element	corotTruss	156	96	36	$A	$matTag
element	corotTruss	157	97	37	$A	$matTag
element	corotTruss	158	98	38	$A	$matTag
element	corotTruss	159	99	39	$A	$matTag
element	corotTruss	160	100	40	$A	$matTag
element	corotTruss	161	101	41	$A	$matTag
element	corotTruss	162	102	42	$A	$matTag
element	corotTruss	163	103	43	$A	$matTag
element	corotTruss	164	104	44	$A	$matTag
element	corotTruss	165	105	45	$A	$matTag
element	corotTruss	166	106	46	$A	$matTag
element	corotTruss	167	107	47	$A	$matTag
element	corotTruss	168	108	48	$A	$matTag
element	corotTruss	169	109	49	$A	$matTag
element	corotTruss	170	110	50	$A	$matTag
element	corotTruss	171	111	51	$A	$matTag
element	corotTruss	172	112	52	$A	$matTag
element	corotTruss	173	113	53	$A	$matTag
element	corotTruss	174	114	54	$A	$matTag
element	corotTruss	175	115	55	$A	$matTag
element	corotTruss	176	116	56	$A	$matTag
element	corotTruss	177	117	57	$A	$matTag
element	corotTruss	178	118	58	$A	$matTag
element	corotTruss	179	119	59	$A	$matTag
element	corotTruss	180	61	30	$A	$matTag
element	corotTruss	181	62	31	$A	$matTag
element	corotTruss	182	63	32	$A	$matTag
element	corotTruss	183	64	33	$A	$matTag
element	corotTruss	184	65	34	$A	$matTag
element	corotTruss	185	66	35	$A	$matTag
element	corotTruss	186	67	36	$A	$matTag
element	corotTruss	187	68	37	$A	$matTag
element	corotTruss	188	69	38	$A	$matTag
element	corotTruss	189	70	39	$A	$matTag
element	corotTruss	190	71	40	$A	$matTag
element	corotTruss	191	72	41	$A	$matTag
element	corotTruss	192	73	42	$A	$matTag
element	corotTruss	193	74	43	$A	$matTag
element	corotTruss	194	75	44	$A	$matTag
element	corotTruss	195	76	45	$A	$matTag
element	corotTruss	196	77	46	$A	$matTag
element	corotTruss	197	78	47	$A	$matTag
element	corotTruss	198	79	48	$A	$matTag
element	corotTruss	199	80	49	$A	$matTag
element	corotTruss	200	81	50	$A	$matTag
element	corotTruss	201	82	51	$A	$matTag
element	corotTruss	202	83	52	$A	$matTag
element	corotTruss	203	84	53	$A	$matTag
element	corotTruss	204	85	54	$A	$matTag
element	corotTruss	205	86	55	$A	$matTag
element	corotTruss	206	87	56	$A	$matTag
element	corotTruss	207	88	57	$A	$matTag
element	corotTruss	208	89	58	$A	$matTag
element	corotTruss	209	60	59	$A	$matTag
element	corotTruss	210	90	60	$A	$matTag
element	corotTruss	211	91	61	$A	$matTag
element	corotTruss	212	92	62	$A	$matTag
element	corotTruss	213	93	63	$A	$matTag
element	corotTruss	214	94	64	$A	$matTag
element	corotTruss	215	95	65	$A	$matTag
element	corotTruss	216	96	66	$A	$matTag
element	corotTruss	217	97	67	$A	$matTag
element	corotTruss	218	98	68	$A	$matTag
element	corotTruss	219	99	69	$A	$matTag
element	corotTruss	220	100	70	$A	$matTag
element	corotTruss	221	101	71	$A	$matTag
element	corotTruss	222	102	72	$A	$matTag
element	corotTruss	223	103	73	$A	$matTag
element	corotTruss	224	104	74	$A	$matTag
element	corotTruss	225	105	75	$A	$matTag
element	corotTruss	226	106	76	$A	$matTag
element	corotTruss	227	107	77	$A	$matTag
element	corotTruss	228	108	78	$A	$matTag
element	corotTruss	229	109	79	$A	$matTag
element	corotTruss	230	110	80	$A	$matTag
element	corotTruss	231	111	81	$A	$matTag
element	corotTruss	232	112	82	$A	$matTag
element	corotTruss	233	113	83	$A	$matTag
element	corotTruss	234	114	84	$A	$matTag
element	corotTruss	235	115	85	$A	$matTag
element	corotTruss	236	116	86	$A	$matTag
element	corotTruss	237	117	87	$A	$matTag
element	corotTruss	238	118	88	$A	$matTag
element	corotTruss	239	119	89	$A	$matTag
element	corotTruss	240	61	60	$A	$matTag
element	corotTruss	241	62	61	$A	$matTag
element	corotTruss	242	63	62	$A	$matTag
element	corotTruss	243	64	63	$A	$matTag
element	corotTruss	244	65	64	$A	$matTag
element	corotTruss	245	66	65	$A	$matTag
element	corotTruss	246	67	66	$A	$matTag
element	corotTruss	247	68	67	$A	$matTag
element	corotTruss	248	69	68	$A	$matTag
element	corotTruss	249	70	69	$A	$matTag
element	corotTruss	250	71	70	$A	$matTag
element	corotTruss	251	72	71	$A	$matTag
element	corotTruss	252	73	72	$A	$matTag
element	corotTruss	253	74	73	$A	$matTag
element	corotTruss	254	75	74	$A	$matTag
element	corotTruss	255	76	75	$A	$matTag
element	corotTruss	256	77	76	$A	$matTag
element	corotTruss	257	78	77	$A	$matTag
element	corotTruss	258	79	78	$A	$matTag
element	corotTruss	259	80	79	$A	$matTag
element	corotTruss	260	81	80	$A	$matTag
element	corotTruss	261	82	81	$A	$matTag
element	corotTruss	262	83	82	$A	$matTag
element	corotTruss	263	84	83	$A	$matTag
element	corotTruss	264	85	84	$A	$matTag
element	corotTruss	265	86	85	$A	$matTag
element	corotTruss	266	87	86	$A	$matTag
element	corotTruss	267	88	87	$A	$matTag
element	corotTruss	268	89	88	$A	$matTag
element	corotTruss	269	60	89	$A	$matTag
element	corotTruss	270	61	90	$A	$matTag
element	corotTruss	271	62	91	$A	$matTag
element	corotTruss	272	63	92	$A	$matTag
element	corotTruss	273	64	93	$A	$matTag
element	corotTruss	274	65	94	$A	$matTag
element	corotTruss	275	66	95	$A	$matTag
element	corotTruss	276	67	96	$A	$matTag
element	corotTruss	277	68	97	$A	$matTag
element	corotTruss	278	69	98	$A	$matTag
element	corotTruss	279	70	99	$A	$matTag
element	corotTruss	280	71	100	$A	$matTag
element	corotTruss	281	72	101	$A	$matTag
element	corotTruss	282	73	102	$A	$matTag
element	corotTruss	283	74	103	$A	$matTag
element	corotTruss	284	75	104	$A	$matTag
element	corotTruss	285	76	105	$A	$matTag
element	corotTruss	286	77	106	$A	$matTag
element	corotTruss	287	78	107	$A	$matTag
element	corotTruss	288	79	108	$A	$matTag
element	corotTruss	289	80	109	$A	$matTag
element	corotTruss	290	81	110	$A	$matTag
element	corotTruss	291	82	111	$A	$matTag
element	corotTruss	292	83	112	$A	$matTag
element	corotTruss	293	84	113	$A	$matTag
element	corotTruss	294	85	114	$A	$matTag
element	corotTruss	295	86	115	$A	$matTag
element	corotTruss	296	87	116	$A	$matTag
element	corotTruss	297	88	117	$A	$matTag
element	corotTruss	298	89	118	$A	$matTag
element	corotTruss	299	60	119	$A	$matTag
element	corotTruss	300	91	90	$A	$matTag
element	corotTruss	301	92	91	$A	$matTag
element	corotTruss	302	93	92	$A	$matTag
element	corotTruss	303	94	93	$A	$matTag
element	corotTruss	304	95	94	$A	$matTag
element	corotTruss	305	96	95	$A	$matTag
element	corotTruss	306	97	96	$A	$matTag
element	corotTruss	307	98	97	$A	$matTag
element	corotTruss	308	99	98	$A	$matTag
element	corotTruss	309	100	99	$A	$matTag
element	corotTruss	310	101	100	$A	$matTag
element	corotTruss	311	102	101	$A	$matTag
element	corotTruss	312	103	102	$A	$matTag
element	corotTruss	313	104	103	$A	$matTag
element	corotTruss	314	105	104	$A	$matTag
element	corotTruss	315	106	105	$A	$matTag
element	corotTruss	316	107	106	$A	$matTag
element	corotTruss	317	108	107	$A	$matTag
element	corotTruss	318	109	108	$A	$matTag
element	corotTruss	319	110	109	$A	$matTag
element	corotTruss	320	111	110	$A	$matTag
element	corotTruss	321	112	111	$A	$matTag
element	corotTruss	322	113	112	$A	$matTag
element	corotTruss	323	114	113	$A	$matTag
element	corotTruss	324	115	114	$A	$matTag
element	corotTruss	325	116	115	$A	$matTag
element	corotTruss	326	117	116	$A	$matTag
element	corotTruss	327	118	117	$A	$matTag
element	corotTruss	328	119	118	$A	$matTag
element	corotTruss	329	90	119	$A	$matTag
element	corotTruss	330	120	90	$A	$matTag
element	corotTruss	331	121	91	$A	$matTag
element	corotTruss	332	122	92	$A	$matTag
element	corotTruss	333	123	93	$A	$matTag
element	corotTruss	334	124	94	$A	$matTag
element	corotTruss	335	125	95	$A	$matTag
element	corotTruss	336	126	96	$A	$matTag
element	corotTruss	337	127	97	$A	$matTag
element	corotTruss	338	128	98	$A	$matTag
element	corotTruss	339	129	99	$A	$matTag
element	corotTruss	340	130	100	$A	$matTag
element	corotTruss	341	131	101	$A	$matTag
element	corotTruss	342	132	102	$A	$matTag
element	corotTruss	343	133	103	$A	$matTag
element	corotTruss	344	134	104	$A	$matTag
element	corotTruss	345	135	105	$A	$matTag
element	corotTruss	346	136	106	$A	$matTag
element	corotTruss	347	137	107	$A	$matTag
element	corotTruss	348	138	108	$A	$matTag
element	corotTruss	349	139	109	$A	$matTag
element	corotTruss	350	140	110	$A	$matTag
element	corotTruss	351	141	111	$A	$matTag
element	corotTruss	352	142	112	$A	$matTag
element	corotTruss	353	143	113	$A	$matTag
element	corotTruss	354	144	114	$A	$matTag
element	corotTruss	355	145	115	$A	$matTag
element	corotTruss	356	146	116	$A	$matTag
element	corotTruss	357	147	117	$A	$matTag
element	corotTruss	358	148	118	$A	$matTag
element	corotTruss	359	149	119	$A	$matTag
element	corotTruss	360	121	90	$A	$matTag
element	corotTruss	361	122	91	$A	$matTag
element	corotTruss	362	123	92	$A	$matTag
element	corotTruss	363	124	93	$A	$matTag
element	corotTruss	364	125	94	$A	$matTag
element	corotTruss	365	126	95	$A	$matTag
element	corotTruss	366	127	96	$A	$matTag
element	corotTruss	367	128	97	$A	$matTag
element	corotTruss	368	129	98	$A	$matTag
element	corotTruss	369	130	99	$A	$matTag
element	corotTruss	370	131	100	$A	$matTag
element	corotTruss	371	132	101	$A	$matTag
element	corotTruss	372	133	102	$A	$matTag
element	corotTruss	373	134	103	$A	$matTag
element	corotTruss	374	135	104	$A	$matTag
element	corotTruss	375	136	105	$A	$matTag
element	corotTruss	376	137	106	$A	$matTag
element	corotTruss	377	138	107	$A	$matTag
element	corotTruss	378	139	108	$A	$matTag
element	corotTruss	379	140	109	$A	$matTag
element	corotTruss	380	141	110	$A	$matTag
element	corotTruss	381	142	111	$A	$matTag
element	corotTruss	382	143	112	$A	$matTag
element	corotTruss	383	144	113	$A	$matTag
element	corotTruss	384	145	114	$A	$matTag
element	corotTruss	385	146	115	$A	$matTag
element	corotTruss	386	147	116	$A	$matTag
element	corotTruss	387	148	117	$A	$matTag
element	corotTruss	388	149	118	$A	$matTag
element	corotTruss	389	120	119	$A	$matTag
element	corotTruss	390	150	90	$A	$matTag
element	corotTruss	391	151	91	$A	$matTag
element	corotTruss	392	152	92	$A	$matTag
element	corotTruss	393	153	93	$A	$matTag
element	corotTruss	394	154	94	$A	$matTag
element	corotTruss	395	155	95	$A	$matTag
element	corotTruss	396	156	96	$A	$matTag
element	corotTruss	397	157	97	$A	$matTag
element	corotTruss	398	158	98	$A	$matTag
element	corotTruss	399	159	99	$A	$matTag
element	corotTruss	400	160	100	$A	$matTag
element	corotTruss	401	161	101	$A	$matTag
element	corotTruss	402	162	102	$A	$matTag
element	corotTruss	403	163	103	$A	$matTag
element	corotTruss	404	164	104	$A	$matTag
element	corotTruss	405	165	105	$A	$matTag
element	corotTruss	406	166	106	$A	$matTag
element	corotTruss	407	167	107	$A	$matTag
element	corotTruss	408	168	108	$A	$matTag
element	corotTruss	409	169	109	$A	$matTag
element	corotTruss	410	170	110	$A	$matTag
element	corotTruss	411	171	111	$A	$matTag
element	corotTruss	412	172	112	$A	$matTag
element	corotTruss	413	173	113	$A	$matTag
element	corotTruss	414	174	114	$A	$matTag
element	corotTruss	415	175	115	$A	$matTag
element	corotTruss	416	176	116	$A	$matTag
element	corotTruss	417	177	117	$A	$matTag
element	corotTruss	418	178	118	$A	$matTag
element	corotTruss	419	179	119	$A	$matTag
element	corotTruss	420	120	60	$A	$matTag
element	corotTruss	421	121	61	$A	$matTag
element	corotTruss	422	122	62	$A	$matTag
element	corotTruss	423	123	63	$A	$matTag
element	corotTruss	424	124	64	$A	$matTag
element	corotTruss	425	125	65	$A	$matTag
element	corotTruss	426	126	66	$A	$matTag
element	corotTruss	427	127	67	$A	$matTag
element	corotTruss	428	128	68	$A	$matTag
element	corotTruss	429	129	69	$A	$matTag
element	corotTruss	430	130	70	$A	$matTag
element	corotTruss	431	131	71	$A	$matTag
element	corotTruss	432	132	72	$A	$matTag
element	corotTruss	433	133	73	$A	$matTag
element	corotTruss	434	134	74	$A	$matTag
element	corotTruss	435	135	75	$A	$matTag
element	corotTruss	436	136	76	$A	$matTag
element	corotTruss	437	137	77	$A	$matTag
element	corotTruss	438	138	78	$A	$matTag
element	corotTruss	439	139	79	$A	$matTag
element	corotTruss	440	140	80	$A	$matTag
element	corotTruss	441	141	81	$A	$matTag
element	corotTruss	442	142	82	$A	$matTag
element	corotTruss	443	143	83	$A	$matTag
element	corotTruss	444	144	84	$A	$matTag
element	corotTruss	445	145	85	$A	$matTag
element	corotTruss	446	146	86	$A	$matTag
element	corotTruss	447	147	87	$A	$matTag
element	corotTruss	448	148	88	$A	$matTag
element	corotTruss	449	149	89	$A	$matTag
element	corotTruss	450	121	120	$A	$matTag
element	corotTruss	451	122	121	$A	$matTag
element	corotTruss	452	123	122	$A	$matTag
element	corotTruss	453	124	123	$A	$matTag
element	corotTruss	454	125	124	$A	$matTag
element	corotTruss	455	126	125	$A	$matTag
element	corotTruss	456	127	126	$A	$matTag
element	corotTruss	457	128	127	$A	$matTag
element	corotTruss	458	129	128	$A	$matTag
element	corotTruss	459	130	129	$A	$matTag
element	corotTruss	460	131	130	$A	$matTag
element	corotTruss	461	132	131	$A	$matTag
element	corotTruss	462	133	132	$A	$matTag
element	corotTruss	463	134	133	$A	$matTag
element	corotTruss	464	135	134	$A	$matTag
element	corotTruss	465	136	135	$A	$matTag
element	corotTruss	466	137	136	$A	$matTag
element	corotTruss	467	138	137	$A	$matTag
element	corotTruss	468	139	138	$A	$matTag
element	corotTruss	469	140	139	$A	$matTag
element	corotTruss	470	141	140	$A	$matTag
element	corotTruss	471	142	141	$A	$matTag
element	corotTruss	472	143	142	$A	$matTag
element	corotTruss	473	144	143	$A	$matTag
element	corotTruss	474	145	144	$A	$matTag
element	corotTruss	475	146	145	$A	$matTag
element	corotTruss	476	147	146	$A	$matTag
element	corotTruss	477	148	147	$A	$matTag
element	corotTruss	478	149	148	$A	$matTag
element	corotTruss	479	120	149	$A	$matTag
element	corotTruss	480	150	120	$A	$matTag
element	corotTruss	481	151	121	$A	$matTag
element	corotTruss	482	152	122	$A	$matTag
element	corotTruss	483	153	123	$A	$matTag
element	corotTruss	484	154	124	$A	$matTag
element	corotTruss	485	155	125	$A	$matTag
element	corotTruss	486	156	126	$A	$matTag
element	corotTruss	487	157	127	$A	$matTag
element	corotTruss	488	158	128	$A	$matTag
element	corotTruss	489	159	129	$A	$matTag
element	corotTruss	490	160	130	$A	$matTag
element	corotTruss	491	161	131	$A	$matTag
element	corotTruss	492	162	132	$A	$matTag
element	corotTruss	493	163	133	$A	$matTag
element	corotTruss	494	164	134	$A	$matTag
element	corotTruss	495	165	135	$A	$matTag
element	corotTruss	496	166	136	$A	$matTag
element	corotTruss	497	167	137	$A	$matTag
element	corotTruss	498	168	138	$A	$matTag
element	corotTruss	499	169	139	$A	$matTag
element	corotTruss	500	170	140	$A	$matTag
element	corotTruss	501	171	141	$A	$matTag
element	corotTruss	502	172	142	$A	$matTag
element	corotTruss	503	173	143	$A	$matTag
element	corotTruss	504	174	144	$A	$matTag
element	corotTruss	505	175	145	$A	$matTag
element	corotTruss	506	176	146	$A	$matTag
element	corotTruss	507	177	147	$A	$matTag
element	corotTruss	508	178	148	$A	$matTag
element	corotTruss	509	179	149	$A	$matTag
element	corotTruss	510	121	150	$A	$matTag
element	corotTruss	511	122	151	$A	$matTag
element	corotTruss	512	123	152	$A	$matTag
element	corotTruss	513	124	153	$A	$matTag
element	corotTruss	514	125	154	$A	$matTag
element	corotTruss	515	126	155	$A	$matTag
element	corotTruss	516	127	156	$A	$matTag
element	corotTruss	517	128	157	$A	$matTag
element	corotTruss	518	129	158	$A	$matTag
element	corotTruss	519	130	159	$A	$matTag
element	corotTruss	520	131	160	$A	$matTag
element	corotTruss	521	132	161	$A	$matTag
element	corotTruss	522	133	162	$A	$matTag
element	corotTruss	523	134	163	$A	$matTag
element	corotTruss	524	135	164	$A	$matTag
element	corotTruss	525	136	165	$A	$matTag
element	corotTruss	526	137	166	$A	$matTag
element	corotTruss	527	138	167	$A	$matTag
element	corotTruss	528	139	168	$A	$matTag
element	corotTruss	529	140	169	$A	$matTag
element	corotTruss	530	141	170	$A	$matTag
element	corotTruss	531	142	171	$A	$matTag
element	corotTruss	532	143	172	$A	$matTag
element	corotTruss	533	144	173	$A	$matTag
element	corotTruss	534	145	174	$A	$matTag
element	corotTruss	535	146	175	$A	$matTag
element	corotTruss	536	147	176	$A	$matTag
element	corotTruss	537	148	177	$A	$matTag
element	corotTruss	538	149	178	$A	$matTag
element	corotTruss	539	120	179	$A	$matTag
element	corotTruss	540	151	150	$A	$matTag
element	corotTruss	541	152	151	$A	$matTag
element	corotTruss	542	153	152	$A	$matTag
element	corotTruss	543	154	153	$A	$matTag
element	corotTruss	544	155	154	$A	$matTag
element	corotTruss	545	156	155	$A	$matTag
element	corotTruss	546	157	156	$A	$matTag
element	corotTruss	547	158	157	$A	$matTag
element	corotTruss	548	159	158	$A	$matTag
element	corotTruss	549	160	159	$A	$matTag
element	corotTruss	550	161	160	$A	$matTag
element	corotTruss	551	162	161	$A	$matTag
element	corotTruss	552	163	162	$A	$matTag
element	corotTruss	553	164	163	$A	$matTag
element	corotTruss	554	165	164	$A	$matTag
element	corotTruss	555	166	165	$A	$matTag
element	corotTruss	556	167	166	$A	$matTag
element	corotTruss	557	168	167	$A	$matTag
element	corotTruss	558	169	168	$A	$matTag
element	corotTruss	559	170	169	$A	$matTag
element	corotTruss	560	171	170	$A	$matTag
element	corotTruss	561	172	171	$A	$matTag
element	corotTruss	562	173	172	$A	$matTag
element	corotTruss	563	174	173	$A	$matTag
element	corotTruss	564	175	174	$A	$matTag
element	corotTruss	565	176	175	$A	$matTag
element	corotTruss	566	177	176	$A	$matTag
element	corotTruss	567	178	177	$A	$matTag
element	corotTruss	568	179	178	$A	$matTag
element	corotTruss	569	150	179	$A	$matTag
element	corotTruss	570	180	150	$A	$matTag
element	corotTruss	571	181	151	$A	$matTag
element	corotTruss	572	182	152	$A	$matTag
element	corotTruss	573	183	153	$A	$matTag
element	corotTruss	574	184	154	$A	$matTag
element	corotTruss	575	185	155	$A	$matTag
element	corotTruss	576	186	156	$A	$matTag
element	corotTruss	577	187	157	$A	$matTag
element	corotTruss	578	188	158	$A	$matTag
element	corotTruss	579	189	159	$A	$matTag
element	corotTruss	580	190	160	$A	$matTag
element	corotTruss	581	191	161	$A	$matTag
element	corotTruss	582	192	162	$A	$matTag
element	corotTruss	583	193	163	$A	$matTag
element	corotTruss	584	194	164	$A	$matTag
element	corotTruss	585	195	165	$A	$matTag
element	corotTruss	586	196	166	$A	$matTag
element	corotTruss	587	197	167	$A	$matTag
element	corotTruss	588	198	168	$A	$matTag
element	corotTruss	589	199	169	$A	$matTag
element	corotTruss	590	200	170	$A	$matTag
element	corotTruss	591	201	171	$A	$matTag
element	corotTruss	592	202	172	$A	$matTag
element	corotTruss	593	203	173	$A	$matTag
element	corotTruss	594	204	174	$A	$matTag
element	corotTruss	595	205	175	$A	$matTag
element	corotTruss	596	206	176	$A	$matTag
element	corotTruss	597	207	177	$A	$matTag
element	corotTruss	598	208	178	$A	$matTag
element	corotTruss	599	209	179	$A	$matTag
element	corotTruss	600	181	150	$A	$matTag
element	corotTruss	601	182	151	$A	$matTag
element	corotTruss	602	183	152	$A	$matTag
element	corotTruss	603	184	153	$A	$matTag
element	corotTruss	604	185	154	$A	$matTag
element	corotTruss	605	186	155	$A	$matTag
element	corotTruss	606	187	156	$A	$matTag
element	corotTruss	607	188	157	$A	$matTag
element	corotTruss	608	189	158	$A	$matTag
element	corotTruss	609	190	159	$A	$matTag
element	corotTruss	610	191	160	$A	$matTag
element	corotTruss	611	192	161	$A	$matTag
element	corotTruss	612	193	162	$A	$matTag
element	corotTruss	613	194	163	$A	$matTag
element	corotTruss	614	195	164	$A	$matTag
element	corotTruss	615	196	165	$A	$matTag
element	corotTruss	616	197	166	$A	$matTag
element	corotTruss	617	198	167	$A	$matTag
element	corotTruss	618	199	168	$A	$matTag
element	corotTruss	619	200	169	$A	$matTag
element	corotTruss	620	201	170	$A	$matTag
element	corotTruss	621	202	171	$A	$matTag
element	corotTruss	622	203	172	$A	$matTag
element	corotTruss	623	204	173	$A	$matTag
element	corotTruss	624	205	174	$A	$matTag
element	corotTruss	625	206	175	$A	$matTag
element	corotTruss	626	207	176	$A	$matTag
element	corotTruss	627	208	177	$A	$matTag
element	corotTruss	628	209	178	$A	$matTag
element	corotTruss	629	180	179	$A	$matTag
element	corotTruss	630	210	150	$A	$matTag
element	corotTruss	631	211	151	$A	$matTag
element	corotTruss	632	212	152	$A	$matTag
element	corotTruss	633	213	153	$A	$matTag
element	corotTruss	634	214	154	$A	$matTag
element	corotTruss	635	215	155	$A	$matTag
element	corotTruss	636	216	156	$A	$matTag
element	corotTruss	637	217	157	$A	$matTag
element	corotTruss	638	218	158	$A	$matTag
element	corotTruss	639	219	159	$A	$matTag
element	corotTruss	640	220	160	$A	$matTag
element	corotTruss	641	221	161	$A	$matTag
element	corotTruss	642	222	162	$A	$matTag
element	corotTruss	643	223	163	$A	$matTag
element	corotTruss	644	224	164	$A	$matTag
element	corotTruss	645	225	165	$A	$matTag
element	corotTruss	646	226	166	$A	$matTag
element	corotTruss	647	227	167	$A	$matTag
element	corotTruss	648	228	168	$A	$matTag
element	corotTruss	649	229	169	$A	$matTag
element	corotTruss	650	230	170	$A	$matTag
element	corotTruss	651	231	171	$A	$matTag
element	corotTruss	652	232	172	$A	$matTag
element	corotTruss	653	233	173	$A	$matTag
element	corotTruss	654	234	174	$A	$matTag
element	corotTruss	655	235	175	$A	$matTag
element	corotTruss	656	236	176	$A	$matTag
element	corotTruss	657	237	177	$A	$matTag
element	corotTruss	658	238	178	$A	$matTag
element	corotTruss	659	239	179	$A	$matTag
element	corotTruss	660	180	120	$A	$matTag
element	corotTruss	661	181	121	$A	$matTag
element	corotTruss	662	182	122	$A	$matTag
element	corotTruss	663	183	123	$A	$matTag
element	corotTruss	664	184	124	$A	$matTag
element	corotTruss	665	185	125	$A	$matTag
element	corotTruss	666	186	126	$A	$matTag
element	corotTruss	667	187	127	$A	$matTag
element	corotTruss	668	188	128	$A	$matTag
element	corotTruss	669	189	129	$A	$matTag
element	corotTruss	670	190	130	$A	$matTag
element	corotTruss	671	191	131	$A	$matTag
element	corotTruss	672	192	132	$A	$matTag
element	corotTruss	673	193	133	$A	$matTag
element	corotTruss	674	194	134	$A	$matTag
element	corotTruss	675	195	135	$A	$matTag
element	corotTruss	676	196	136	$A	$matTag
element	corotTruss	677	197	137	$A	$matTag
element	corotTruss	678	198	138	$A	$matTag
element	corotTruss	679	199	139	$A	$matTag
element	corotTruss	680	200	140	$A	$matTag
element	corotTruss	681	201	141	$A	$matTag
element	corotTruss	682	202	142	$A	$matTag
element	corotTruss	683	203	143	$A	$matTag
element	corotTruss	684	204	144	$A	$matTag
element	corotTruss	685	205	145	$A	$matTag
element	corotTruss	686	206	146	$A	$matTag
element	corotTruss	687	207	147	$A	$matTag
element	corotTruss	688	208	148	$A	$matTag
element	corotTruss	689	209	149	$A	$matTag
element	corotTruss	690	181	180	$A	$matTag
element	corotTruss	691	182	181	$A	$matTag
element	corotTruss	692	183	182	$A	$matTag
element	corotTruss	693	184	183	$A	$matTag
element	corotTruss	694	185	184	$A	$matTag
element	corotTruss	695	186	185	$A	$matTag
element	corotTruss	696	187	186	$A	$matTag
element	corotTruss	697	188	187	$A	$matTag
element	corotTruss	698	189	188	$A	$matTag
element	corotTruss	699	190	189	$A	$matTag
element	corotTruss	700	191	190	$A	$matTag
element	corotTruss	701	192	191	$A	$matTag
element	corotTruss	702	193	192	$A	$matTag
element	corotTruss	703	194	193	$A	$matTag
element	corotTruss	704	195	194	$A	$matTag
element	corotTruss	705	196	195	$A	$matTag
element	corotTruss	706	197	196	$A	$matTag
element	corotTruss	707	198	197	$A	$matTag
element	corotTruss	708	199	198	$A	$matTag
element	corotTruss	709	200	199	$A	$matTag
element	corotTruss	710	201	200	$A	$matTag
element	corotTruss	711	202	201	$A	$matTag
element	corotTruss	712	203	202	$A	$matTag
element	corotTruss	713	204	203	$A	$matTag
element	corotTruss	714	205	204	$A	$matTag
element	corotTruss	715	206	205	$A	$matTag
element	corotTruss	716	207	206	$A	$matTag
element	corotTruss	717	208	207	$A	$matTag
element	corotTruss	718	209	208	$A	$matTag
element	corotTruss	719	180	209	$A	$matTag
element	corotTruss	720	210	180	$A	$matTag
element	corotTruss	721	211	181	$A	$matTag
element	corotTruss	722	212	182	$A	$matTag
element	corotTruss	723	213	183	$A	$matTag
element	corotTruss	724	214	184	$A	$matTag
element	corotTruss	725	215	185	$A	$matTag
element	corotTruss	726	216	186	$A	$matTag
element	corotTruss	727	217	187	$A	$matTag
element	corotTruss	728	218	188	$A	$matTag
element	corotTruss	729	219	189	$A	$matTag
element	corotTruss	730	220	190	$A	$matTag
element	corotTruss	731	221	191	$A	$matTag
element	corotTruss	732	222	192	$A	$matTag
element	corotTruss	733	223	193	$A	$matTag
element	corotTruss	734	224	194	$A	$matTag
element	corotTruss	735	225	195	$A	$matTag
element	corotTruss	736	226	196	$A	$matTag
element	corotTruss	737	227	197	$A	$matTag
element	corotTruss	738	228	198	$A	$matTag
element	corotTruss	739	229	199	$A	$matTag
element	corotTruss	740	230	200	$A	$matTag
element	corotTruss	741	231	201	$A	$matTag
element	corotTruss	742	232	202	$A	$matTag
element	corotTruss	743	233	203	$A	$matTag
element	corotTruss	744	234	204	$A	$matTag
element	corotTruss	745	235	205	$A	$matTag
element	corotTruss	746	236	206	$A	$matTag
element	corotTruss	747	237	207	$A	$matTag
element	corotTruss	748	238	208	$A	$matTag
element	corotTruss	749	239	209	$A	$matTag
element	corotTruss	750	181	210	$A	$matTag
element	corotTruss	751	182	211	$A	$matTag
element	corotTruss	752	183	212	$A	$matTag
element	corotTruss	753	184	213	$A	$matTag
element	corotTruss	754	185	214	$A	$matTag
element	corotTruss	755	186	215	$A	$matTag
element	corotTruss	756	187	216	$A	$matTag
element	corotTruss	757	188	217	$A	$matTag
element	corotTruss	758	189	218	$A	$matTag
element	corotTruss	759	190	219	$A	$matTag
element	corotTruss	760	191	220	$A	$matTag
element	corotTruss	761	192	221	$A	$matTag
element	corotTruss	762	193	222	$A	$matTag
element	corotTruss	763	194	223	$A	$matTag
element	corotTruss	764	195	224	$A	$matTag
element	corotTruss	765	196	225	$A	$matTag
element	corotTruss	766	197	226	$A	$matTag
element	corotTruss	767	198	227	$A	$matTag
element	corotTruss	768	199	228	$A	$matTag
element	corotTruss	769	200	229	$A	$matTag
element	corotTruss	770	201	230	$A	$matTag
element	corotTruss	771	202	231	$A	$matTag
element	corotTruss	772	203	232	$A	$matTag
element	corotTruss	773	204	233	$A	$matTag
element	corotTruss	774	205	234	$A	$matTag
element	corotTruss	775	206	235	$A	$matTag
element	corotTruss	776	207	236	$A	$matTag
element	corotTruss	777	208	237	$A	$matTag
element	corotTruss	778	209	238	$A	$matTag
element	corotTruss	779	180	239	$A	$matTag
element	corotTruss	780	211	210	$A	$matTag
element	corotTruss	781	212	211	$A	$matTag
element	corotTruss	782	213	212	$A	$matTag
element	corotTruss	783	214	213	$A	$matTag
element	corotTruss	784	215	214	$A	$matTag
element	corotTruss	785	216	215	$A	$matTag
element	corotTruss	786	217	216	$A	$matTag
element	corotTruss	787	218	217	$A	$matTag
element	corotTruss	788	219	218	$A	$matTag
element	corotTruss	789	220	219	$A	$matTag
element	corotTruss	790	221	220	$A	$matTag
element	corotTruss	791	222	221	$A	$matTag
element	corotTruss	792	223	222	$A	$matTag
element	corotTruss	793	224	223	$A	$matTag
element	corotTruss	794	225	224	$A	$matTag
element	corotTruss	795	226	225	$A	$matTag
element	corotTruss	796	227	226	$A	$matTag
element	corotTruss	797	228	227	$A	$matTag
element	corotTruss	798	229	228	$A	$matTag
element	corotTruss	799	230	229	$A	$matTag
element	corotTruss	800	231	230	$A	$matTag
element	corotTruss	801	232	231	$A	$matTag
element	corotTruss	802	233	232	$A	$matTag
element	corotTruss	803	234	233	$A	$matTag
element	corotTruss	804	235	234	$A	$matTag
element	corotTruss	805	236	235	$A	$matTag
element	corotTruss	806	237	236	$A	$matTag
element	corotTruss	807	238	237	$A	$matTag
element	corotTruss	808	239	238	$A	$matTag
element	corotTruss	809	210	239	$A	$matTag
element	corotTruss	810	240	210	$A	$matTag
element	corotTruss	811	241	211	$A	$matTag
element	corotTruss	812	242	212	$A	$matTag
element	corotTruss	813	243	213	$A	$matTag
element	corotTruss	814	244	214	$A	$matTag
element	corotTruss	815	245	215	$A	$matTag
element	corotTruss	816	246	216	$A	$matTag
element	corotTruss	817	247	217	$A	$matTag
element	corotTruss	818	248	218	$A	$matTag
element	corotTruss	819	249	219	$A	$matTag
element	corotTruss	820	250	220	$A	$matTag
element	corotTruss	821	251	221	$A	$matTag
element	corotTruss	822	252	222	$A	$matTag
element	corotTruss	823	253	223	$A	$matTag
element	corotTruss	824	254	224	$A	$matTag
element	corotTruss	825	255	225	$A	$matTag
element	corotTruss	826	256	226	$A	$matTag
element	corotTruss	827	257	227	$A	$matTag
element	corotTruss	828	258	228	$A	$matTag
element	corotTruss	829	259	229	$A	$matTag
element	corotTruss	830	260	230	$A	$matTag
element	corotTruss	831	261	231	$A	$matTag
element	corotTruss	832	262	232	$A	$matTag
element	corotTruss	833	263	233	$A	$matTag
element	corotTruss	834	264	234	$A	$matTag
element	corotTruss	835	265	235	$A	$matTag
element	corotTruss	836	266	236	$A	$matTag
element	corotTruss	837	267	237	$A	$matTag
element	corotTruss	838	268	238	$A	$matTag
element	corotTruss	839	269	239	$A	$matTag
element	corotTruss	840	241	210	$A	$matTag
element	corotTruss	841	242	211	$A	$matTag
element	corotTruss	842	243	212	$A	$matTag
element	corotTruss	843	244	213	$A	$matTag
element	corotTruss	844	245	214	$A	$matTag
element	corotTruss	845	246	215	$A	$matTag
element	corotTruss	846	247	216	$A	$matTag
element	corotTruss	847	248	217	$A	$matTag
element	corotTruss	848	249	218	$A	$matTag
element	corotTruss	849	250	219	$A	$matTag
element	corotTruss	850	251	220	$A	$matTag
element	corotTruss	851	252	221	$A	$matTag
element	corotTruss	852	253	222	$A	$matTag
element	corotTruss	853	254	223	$A	$matTag
element	corotTruss	854	255	224	$A	$matTag
element	corotTruss	855	256	225	$A	$matTag
element	corotTruss	856	257	226	$A	$matTag
element	corotTruss	857	258	227	$A	$matTag
element	corotTruss	858	259	228	$A	$matTag
element	corotTruss	859	260	229	$A	$matTag
element	corotTruss	860	261	230	$A	$matTag
element	corotTruss	861	262	231	$A	$matTag
element	corotTruss	862	263	232	$A	$matTag
element	corotTruss	863	264	233	$A	$matTag
element	corotTruss	864	265	234	$A	$matTag
element	corotTruss	865	266	235	$A	$matTag
element	corotTruss	866	267	236	$A	$matTag
element	corotTruss	867	268	237	$A	$matTag
element	corotTruss	868	269	238	$A	$matTag
element	corotTruss	869	240	239	$A	$matTag
element	corotTruss	870	240	180	$A	$matTag
element	corotTruss	871	241	181	$A	$matTag
element	corotTruss	872	242	182	$A	$matTag
element	corotTruss	873	243	183	$A	$matTag
element	corotTruss	874	244	184	$A	$matTag
element	corotTruss	875	245	185	$A	$matTag
element	corotTruss	876	246	186	$A	$matTag
element	corotTruss	877	247	187	$A	$matTag
element	corotTruss	878	248	188	$A	$matTag
element	corotTruss	879	249	189	$A	$matTag
element	corotTruss	880	250	190	$A	$matTag
element	corotTruss	881	251	191	$A	$matTag
element	corotTruss	882	252	192	$A	$matTag
element	corotTruss	883	253	193	$A	$matTag
element	corotTruss	884	254	194	$A	$matTag
element	corotTruss	885	255	195	$A	$matTag
element	corotTruss	886	256	196	$A	$matTag
element	corotTruss	887	257	197	$A	$matTag
element	corotTruss	888	258	198	$A	$matTag
element	corotTruss	889	259	199	$A	$matTag
element	corotTruss	890	260	200	$A	$matTag
element	corotTruss	891	261	201	$A	$matTag
element	corotTruss	892	262	202	$A	$matTag
element	corotTruss	893	263	203	$A	$matTag
element	corotTruss	894	264	204	$A	$matTag
element	corotTruss	895	265	205	$A	$matTag
element	corotTruss	896	266	206	$A	$matTag
element	corotTruss	897	267	207	$A	$matTag
element	corotTruss	898	268	208	$A	$matTag
element	corotTruss	899	269	209	$A	$matTag
element	corotTruss	900	270	210	$A	$matTag
element	corotTruss	901	271	211	$A	$matTag
element	corotTruss	902	272	212	$A	$matTag
element	corotTruss	903	273	213	$A	$matTag
element	corotTruss	904	274	214	$A	$matTag
element	corotTruss	905	275	215	$A	$matTag
element	corotTruss	906	276	216	$A	$matTag
element	corotTruss	907	277	217	$A	$matTag
element	corotTruss	908	278	218	$A	$matTag
element	corotTruss	909	279	219	$A	$matTag
element	corotTruss	910	280	220	$A	$matTag
element	corotTruss	911	281	221	$A	$matTag
element	corotTruss	912	282	222	$A	$matTag
element	corotTruss	913	283	223	$A	$matTag
element	corotTruss	914	284	224	$A	$matTag
element	corotTruss	915	285	225	$A	$matTag
element	corotTruss	916	286	226	$A	$matTag
element	corotTruss	917	287	227	$A	$matTag
element	corotTruss	918	288	228	$A	$matTag
element	corotTruss	919	289	229	$A	$matTag
element	corotTruss	920	290	230	$A	$matTag
element	corotTruss	921	291	231	$A	$matTag
element	corotTruss	922	292	232	$A	$matTag
element	corotTruss	923	293	233	$A	$matTag
element	corotTruss	924	294	234	$A	$matTag
element	corotTruss	925	295	235	$A	$matTag
element	corotTruss	926	296	236	$A	$matTag
element	corotTruss	927	297	237	$A	$matTag
element	corotTruss	928	298	238	$A	$matTag
element	corotTruss	929	299	239	$A	$matTag
element	corotTruss	930	240	270	$A	$matTag
element	corotTruss	931	241	271	$A	$matTag
element	corotTruss	932	242	272	$A	$matTag
element	corotTruss	933	243	273	$A	$matTag
element	corotTruss	934	244	274	$A	$matTag
element	corotTruss	935	245	275	$A	$matTag
element	corotTruss	936	246	276	$A	$matTag
element	corotTruss	937	247	277	$A	$matTag
element	corotTruss	938	248	278	$A	$matTag
element	corotTruss	939	249	279	$A	$matTag
element	corotTruss	940	250	280	$A	$matTag
element	corotTruss	941	251	281	$A	$matTag
element	corotTruss	942	252	282	$A	$matTag
element	corotTruss	943	253	283	$A	$matTag
element	corotTruss	944	254	284	$A	$matTag
element	corotTruss	945	255	285	$A	$matTag
element	corotTruss	946	256	286	$A	$matTag
element	corotTruss	947	257	287	$A	$matTag
element	corotTruss	948	258	288	$A	$matTag
element	corotTruss	949	259	289	$A	$matTag
element	corotTruss	950	260	290	$A	$matTag
element	corotTruss	951	261	291	$A	$matTag
element	corotTruss	952	262	292	$A	$matTag
element	corotTruss	953	263	293	$A	$matTag
element	corotTruss	954	264	294	$A	$matTag
element	corotTruss	955	265	295	$A	$matTag
element	corotTruss	956	266	296	$A	$matTag
element	corotTruss	957	267	297	$A	$matTag
element	corotTruss	958	268	298	$A	$matTag
element	corotTruss	959	269	299	$A	$matTag
element	corotTruss	960	241	270	$A	$matTag
element	corotTruss	961	242	271	$A	$matTag
element	corotTruss	962	243	272	$A	$matTag
element	corotTruss	963	244	273	$A	$matTag
element	corotTruss	964	245	274	$A	$matTag
element	corotTruss	965	246	275	$A	$matTag
element	corotTruss	966	247	276	$A	$matTag
element	corotTruss	967	248	277	$A	$matTag
element	corotTruss	968	249	278	$A	$matTag
element	corotTruss	969	250	279	$A	$matTag
element	corotTruss	970	251	280	$A	$matTag
element	corotTruss	971	252	281	$A	$matTag
element	corotTruss	972	253	282	$A	$matTag
element	corotTruss	973	254	283	$A	$matTag
element	corotTruss	974	255	284	$A	$matTag
element	corotTruss	975	256	285	$A	$matTag
element	corotTruss	976	257	286	$A	$matTag
element	corotTruss	977	258	287	$A	$matTag
element	corotTruss	978	259	288	$A	$matTag
element	corotTruss	979	260	289	$A	$matTag
element	corotTruss	980	261	290	$A	$matTag
element	corotTruss	981	262	291	$A	$matTag
element	corotTruss	982	263	292	$A	$matTag
element	corotTruss	983	264	293	$A	$matTag
element	corotTruss	984	265	294	$A	$matTag
element	corotTruss	985	266	295	$A	$matTag
element	corotTruss	986	267	296	$A	$matTag
element	corotTruss	987	268	297	$A	$matTag
element	corotTruss	988	269	298	$A	$matTag
element	corotTruss	989	240	299	$A	$matTag
element	corotTruss	990	241	240	$A	$matTag
element	corotTruss	991	242	241	$A	$matTag
element	corotTruss	992	243	242	$A	$matTag
element	corotTruss	993	244	243	$A	$matTag
element	corotTruss	994	245	244	$A	$matTag
element	corotTruss	995	246	245	$A	$matTag
element	corotTruss	996	247	246	$A	$matTag
element	corotTruss	997	248	247	$A	$matTag
element	corotTruss	998	249	248	$A	$matTag
element	corotTruss	999	250	249	$A	$matTag
element	corotTruss	1000	251	250	$A	$matTag
element	corotTruss	1001	252	251	$A	$matTag
element	corotTruss	1002	253	252	$A	$matTag
element	corotTruss	1003	254	253	$A	$matTag
element	corotTruss	1004	255	254	$A	$matTag
element	corotTruss	1005	256	255	$A	$matTag
element	corotTruss	1006	257	256	$A	$matTag
element	corotTruss	1007	258	257	$A	$matTag
element	corotTruss	1008	259	258	$A	$matTag
element	corotTruss	1009	260	259	$A	$matTag
element	corotTruss	1010	261	260	$A	$matTag
element	corotTruss	1011	262	261	$A	$matTag
element	corotTruss	1012	263	262	$A	$matTag
element	corotTruss	1013	264	263	$A	$matTag
element	corotTruss	1014	265	264	$A	$matTag
element	corotTruss	1015	266	265	$A	$matTag
element	corotTruss	1016	267	266	$A	$matTag
element	corotTruss	1017	268	267	$A	$matTag
element	corotTruss	1018	269	268	$A	$matTag
element	corotTruss	1019	240	269	$A	$matTag
element	corotTruss	1020	271	270	$A	$matTag
element	corotTruss	1021	272	271	$A	$matTag
element	corotTruss	1022	273	272	$A	$matTag
element	corotTruss	1023	274	273	$A	$matTag
element	corotTruss	1024	275	274	$A	$matTag
element	corotTruss	1025	276	275	$A	$matTag
element	corotTruss	1026	277	276	$A	$matTag
element	corotTruss	1027	278	277	$A	$matTag
element	corotTruss	1028	279	278	$A	$matTag
element	corotTruss	1029	280	279	$A	$matTag
element	corotTruss	1030	281	280	$A	$matTag
element	corotTruss	1031	282	281	$A	$matTag
element	corotTruss	1032	283	282	$A	$matTag
element	corotTruss	1033	284	283	$A	$matTag
element	corotTruss	1034	285	284	$A	$matTag
element	corotTruss	1035	286	285	$A	$matTag
element	corotTruss	1036	287	286	$A	$matTag
element	corotTruss	1037	288	287	$A	$matTag
element	corotTruss	1038	289	288	$A	$matTag
element	corotTruss	1039	290	289	$A	$matTag
element	corotTruss	1040	291	290	$A	$matTag
element	corotTruss	1041	292	291	$A	$matTag
element	corotTruss	1042	293	292	$A	$matTag
element	corotTruss	1043	294	293	$A	$matTag
element	corotTruss	1044	295	294	$A	$matTag
element	corotTruss	1045	296	295	$A	$matTag
element	corotTruss	1046	297	296	$A	$matTag
element	corotTruss	1047	298	297	$A	$matTag
element	corotTruss	1048	299	298	$A	$matTag
element	corotTruss	1049	270	299	$A	$matTag
element	corotTruss	1050	300	270	$A	$matTag
element	corotTruss	1051	301	271	$A	$matTag
element	corotTruss	1052	302	272	$A	$matTag
element	corotTruss	1053	303	273	$A	$matTag
element	corotTruss	1054	304	274	$A	$matTag
element	corotTruss	1055	305	275	$A	$matTag
element	corotTruss	1056	306	276	$A	$matTag
element	corotTruss	1057	307	277	$A	$matTag
element	corotTruss	1058	308	278	$A	$matTag
element	corotTruss	1059	309	279	$A	$matTag
element	corotTruss	1060	310	280	$A	$matTag
element	corotTruss	1061	311	281	$A	$matTag
element	corotTruss	1062	312	282	$A	$matTag
element	corotTruss	1063	313	283	$A	$matTag
element	corotTruss	1064	314	284	$A	$matTag
element	corotTruss	1065	315	285	$A	$matTag
element	corotTruss	1066	316	286	$A	$matTag
element	corotTruss	1067	317	287	$A	$matTag
element	corotTruss	1068	318	288	$A	$matTag
element	corotTruss	1069	319	289	$A	$matTag
element	corotTruss	1070	320	290	$A	$matTag
element	corotTruss	1071	321	291	$A	$matTag
element	corotTruss	1072	322	292	$A	$matTag
element	corotTruss	1073	323	293	$A	$matTag
element	corotTruss	1074	324	294	$A	$matTag
element	corotTruss	1075	325	295	$A	$matTag
element	corotTruss	1076	326	296	$A	$matTag
element	corotTruss	1077	327	297	$A	$matTag
element	corotTruss	1078	328	298	$A	$matTag
element	corotTruss	1079	329	299	$A	$matTag
element	corotTruss	1080	301	270	$A	$matTag
element	corotTruss	1081	302	271	$A	$matTag
element	corotTruss	1082	303	272	$A	$matTag
element	corotTruss	1083	304	273	$A	$matTag
element	corotTruss	1084	305	274	$A	$matTag
element	corotTruss	1085	306	275	$A	$matTag
element	corotTruss	1086	307	276	$A	$matTag
element	corotTruss	1087	308	277	$A	$matTag
element	corotTruss	1088	309	278	$A	$matTag
element	corotTruss	1089	310	279	$A	$matTag
element	corotTruss	1090	311	280	$A	$matTag
element	corotTruss	1091	312	281	$A	$matTag
element	corotTruss	1092	313	282	$A	$matTag
element	corotTruss	1093	314	283	$A	$matTag
element	corotTruss	1094	315	284	$A	$matTag
element	corotTruss	1095	316	285	$A	$matTag
element	corotTruss	1096	317	286	$A	$matTag
element	corotTruss	1097	318	287	$A	$matTag
element	corotTruss	1098	319	288	$A	$matTag
element	corotTruss	1099	320	289	$A	$matTag
element	corotTruss	1100	321	290	$A	$matTag
element	corotTruss	1101	322	291	$A	$matTag
element	corotTruss	1102	323	292	$A	$matTag
element	corotTruss	1103	324	293	$A	$matTag
element	corotTruss	1104	325	294	$A	$matTag
element	corotTruss	1105	326	295	$A	$matTag
element	corotTruss	1106	327	296	$A	$matTag
element	corotTruss	1107	328	297	$A	$matTag
element	corotTruss	1108	329	298	$A	$matTag
element	corotTruss	1109	300	299	$A	$matTag
element	corotTruss	1110	300	240	$A	$matTag
element	corotTruss	1111	301	241	$A	$matTag
element	corotTruss	1112	302	242	$A	$matTag
element	corotTruss	1113	303	243	$A	$matTag
element	corotTruss	1114	304	244	$A	$matTag
element	corotTruss	1115	305	245	$A	$matTag
element	corotTruss	1116	306	246	$A	$matTag
element	corotTruss	1117	307	247	$A	$matTag
element	corotTruss	1118	308	248	$A	$matTag
element	corotTruss	1119	309	249	$A	$matTag
element	corotTruss	1120	310	250	$A	$matTag
element	corotTruss	1121	311	251	$A	$matTag
element	corotTruss	1122	312	252	$A	$matTag
element	corotTruss	1123	313	253	$A	$matTag
element	corotTruss	1124	314	254	$A	$matTag
element	corotTruss	1125	315	255	$A	$matTag
element	corotTruss	1126	316	256	$A	$matTag
element	corotTruss	1127	317	257	$A	$matTag
element	corotTruss	1128	318	258	$A	$matTag
element	corotTruss	1129	319	259	$A	$matTag
element	corotTruss	1130	320	260	$A	$matTag
element	corotTruss	1131	321	261	$A	$matTag
element	corotTruss	1132	322	262	$A	$matTag
element	corotTruss	1133	323	263	$A	$matTag
element	corotTruss	1134	324	264	$A	$matTag
element	corotTruss	1135	325	265	$A	$matTag
element	corotTruss	1136	326	266	$A	$matTag
element	corotTruss	1137	327	267	$A	$matTag
element	corotTruss	1138	328	268	$A	$matTag
element	corotTruss	1139	329	269	$A	$matTag
element	corotTruss	1140	330	270	$A	$matTag
element	corotTruss	1141	331	271	$A	$matTag
element	corotTruss	1142	332	272	$A	$matTag
element	corotTruss	1143	333	273	$A	$matTag
element	corotTruss	1144	334	274	$A	$matTag
element	corotTruss	1145	335	275	$A	$matTag
element	corotTruss	1146	336	276	$A	$matTag
element	corotTruss	1147	337	277	$A	$matTag
element	corotTruss	1148	338	278	$A	$matTag
element	corotTruss	1149	339	279	$A	$matTag
element	corotTruss	1150	340	280	$A	$matTag
element	corotTruss	1151	341	281	$A	$matTag
element	corotTruss	1152	342	282	$A	$matTag
element	corotTruss	1153	343	283	$A	$matTag
element	corotTruss	1154	344	284	$A	$matTag
element	corotTruss	1155	345	285	$A	$matTag
element	corotTruss	1156	346	286	$A	$matTag
element	corotTruss	1157	347	287	$A	$matTag
element	corotTruss	1158	348	288	$A	$matTag
element	corotTruss	1159	349	289	$A	$matTag
element	corotTruss	1160	350	290	$A	$matTag
element	corotTruss	1161	351	291	$A	$matTag
element	corotTruss	1162	352	292	$A	$matTag
element	corotTruss	1163	353	293	$A	$matTag
element	corotTruss	1164	354	294	$A	$matTag
element	corotTruss	1165	355	295	$A	$matTag
element	corotTruss	1166	356	296	$A	$matTag
element	corotTruss	1167	357	297	$A	$matTag
element	corotTruss	1168	358	298	$A	$matTag
element	corotTruss	1169	359	299	$A	$matTag
element	corotTruss	1170	330	300	$A	$matTag
element	corotTruss	1171	331	301	$A	$matTag
element	corotTruss	1172	332	302	$A	$matTag
element	corotTruss	1173	333	303	$A	$matTag
element	corotTruss	1174	334	304	$A	$matTag
element	corotTruss	1175	335	305	$A	$matTag
element	corotTruss	1176	336	306	$A	$matTag
element	corotTruss	1177	337	307	$A	$matTag
element	corotTruss	1178	338	308	$A	$matTag
element	corotTruss	1179	339	309	$A	$matTag
element	corotTruss	1180	340	310	$A	$matTag
element	corotTruss	1181	341	311	$A	$matTag
element	corotTruss	1182	342	312	$A	$matTag
element	corotTruss	1183	343	313	$A	$matTag
element	corotTruss	1184	344	314	$A	$matTag
element	corotTruss	1185	345	315	$A	$matTag
element	corotTruss	1186	346	316	$A	$matTag
element	corotTruss	1187	347	317	$A	$matTag
element	corotTruss	1188	348	318	$A	$matTag
element	corotTruss	1189	349	319	$A	$matTag
element	corotTruss	1190	350	320	$A	$matTag
element	corotTruss	1191	351	321	$A	$matTag
element	corotTruss	1192	352	322	$A	$matTag
element	corotTruss	1193	353	323	$A	$matTag
element	corotTruss	1194	354	324	$A	$matTag
element	corotTruss	1195	355	325	$A	$matTag
element	corotTruss	1196	356	326	$A	$matTag
element	corotTruss	1197	357	327	$A	$matTag
element	corotTruss	1198	358	328	$A	$matTag
element	corotTruss	1199	359	329	$A	$matTag
element	corotTruss	1200	330	301	$A	$matTag
element	corotTruss	1201	331	302	$A	$matTag
element	corotTruss	1202	332	303	$A	$matTag
element	corotTruss	1203	333	304	$A	$matTag
element	corotTruss	1204	334	305	$A	$matTag
element	corotTruss	1205	335	306	$A	$matTag
element	corotTruss	1206	336	307	$A	$matTag
element	corotTruss	1207	337	308	$A	$matTag
element	corotTruss	1208	338	309	$A	$matTag
element	corotTruss	1209	339	310	$A	$matTag
element	corotTruss	1210	340	311	$A	$matTag
element	corotTruss	1211	341	312	$A	$matTag
element	corotTruss	1212	342	313	$A	$matTag
element	corotTruss	1213	343	314	$A	$matTag
element	corotTruss	1214	344	315	$A	$matTag
element	corotTruss	1215	345	316	$A	$matTag
element	corotTruss	1216	346	317	$A	$matTag
element	corotTruss	1217	347	318	$A	$matTag
element	corotTruss	1218	348	319	$A	$matTag
element	corotTruss	1219	349	320	$A	$matTag
element	corotTruss	1220	350	321	$A	$matTag
element	corotTruss	1221	351	322	$A	$matTag
element	corotTruss	1222	352	323	$A	$matTag
element	corotTruss	1223	353	324	$A	$matTag
element	corotTruss	1224	354	325	$A	$matTag
element	corotTruss	1225	355	326	$A	$matTag
element	corotTruss	1226	356	327	$A	$matTag
element	corotTruss	1227	357	328	$A	$matTag
element	corotTruss	1228	358	329	$A	$matTag
element	corotTruss	1229	359	300	$A	$matTag
element	corotTruss	1230	301	300	$A	$matTag
element	corotTruss	1231	302	301	$A	$matTag
element	corotTruss	1232	303	302	$A	$matTag
element	corotTruss	1233	304	303	$A	$matTag
element	corotTruss	1234	305	304	$A	$matTag
element	corotTruss	1235	306	305	$A	$matTag
element	corotTruss	1236	307	306	$A	$matTag
element	corotTruss	1237	308	307	$A	$matTag
element	corotTruss	1238	309	308	$A	$matTag
element	corotTruss	1239	310	309	$A	$matTag
element	corotTruss	1240	311	310	$A	$matTag
element	corotTruss	1241	312	311	$A	$matTag
element	corotTruss	1242	313	312	$A	$matTag
element	corotTruss	1243	314	313	$A	$matTag
element	corotTruss	1244	315	314	$A	$matTag
element	corotTruss	1245	316	315	$A	$matTag
element	corotTruss	1246	317	316	$A	$matTag
element	corotTruss	1247	318	317	$A	$matTag
element	corotTruss	1248	319	318	$A	$matTag
element	corotTruss	1249	320	319	$A	$matTag
element	corotTruss	1250	321	320	$A	$matTag
element	corotTruss	1251	322	321	$A	$matTag
element	corotTruss	1252	323	322	$A	$matTag
element	corotTruss	1253	324	323	$A	$matTag
element	corotTruss	1254	325	324	$A	$matTag
element	corotTruss	1255	326	325	$A	$matTag
element	corotTruss	1256	327	326	$A	$matTag
element	corotTruss	1257	328	327	$A	$matTag
element	corotTruss	1258	329	328	$A	$matTag
element	corotTruss	1259	300	329	$A	$matTag
element	corotTruss	1260	361	330	$A	$matTag
element	corotTruss	1261	362	331	$A	$matTag
element	corotTruss	1262	363	332	$A	$matTag
element	corotTruss	1263	364	333	$A	$matTag
element	corotTruss	1264	365	334	$A	$matTag
element	corotTruss	1265	366	335	$A	$matTag
element	corotTruss	1266	367	336	$A	$matTag
element	corotTruss	1267	368	337	$A	$matTag
element	corotTruss	1268	369	338	$A	$matTag
element	corotTruss	1269	370	339	$A	$matTag
element	corotTruss	1270	371	340	$A	$matTag
element	corotTruss	1271	372	341	$A	$matTag
element	corotTruss	1272	373	342	$A	$matTag
element	corotTruss	1273	374	343	$A	$matTag
element	corotTruss	1274	375	344	$A	$matTag
element	corotTruss	1275	376	345	$A	$matTag
element	corotTruss	1276	377	346	$A	$matTag
element	corotTruss	1277	378	347	$A	$matTag
element	corotTruss	1278	379	348	$A	$matTag
element	corotTruss	1279	380	349	$A	$matTag
element	corotTruss	1280	381	350	$A	$matTag
element	corotTruss	1281	382	351	$A	$matTag
element	corotTruss	1282	383	352	$A	$matTag
element	corotTruss	1283	384	353	$A	$matTag
element	corotTruss	1284	385	354	$A	$matTag
element	corotTruss	1285	386	355	$A	$matTag
element	corotTruss	1286	387	356	$A	$matTag
element	corotTruss	1287	388	357	$A	$matTag
element	corotTruss	1288	389	358	$A	$matTag
element	corotTruss	1289	360	359	$A	$matTag
element	corotTruss	1290	360	330	$A	$matTag
element	corotTruss	1291	361	331	$A	$matTag
element	corotTruss	1292	362	332	$A	$matTag
element	corotTruss	1293	363	333	$A	$matTag
element	corotTruss	1294	364	334	$A	$matTag
element	corotTruss	1295	365	335	$A	$matTag
element	corotTruss	1296	366	336	$A	$matTag
element	corotTruss	1297	367	337	$A	$matTag
element	corotTruss	1298	368	338	$A	$matTag
element	corotTruss	1299	369	339	$A	$matTag
element	corotTruss	1300	370	340	$A	$matTag
element	corotTruss	1301	371	341	$A	$matTag
element	corotTruss	1302	372	342	$A	$matTag
element	corotTruss	1303	373	343	$A	$matTag
element	corotTruss	1304	374	344	$A	$matTag
element	corotTruss	1305	375	345	$A	$matTag
element	corotTruss	1306	376	346	$A	$matTag
element	corotTruss	1307	377	347	$A	$matTag
element	corotTruss	1308	378	348	$A	$matTag
element	corotTruss	1309	379	349	$A	$matTag
element	corotTruss	1310	380	350	$A	$matTag
element	corotTruss	1311	381	351	$A	$matTag
element	corotTruss	1312	382	352	$A	$matTag
element	corotTruss	1313	383	353	$A	$matTag
element	corotTruss	1314	384	354	$A	$matTag
element	corotTruss	1315	385	355	$A	$matTag
element	corotTruss	1316	386	356	$A	$matTag
element	corotTruss	1317	387	357	$A	$matTag
element	corotTruss	1318	388	358	$A	$matTag
element	corotTruss	1319	389	359	$A	$matTag
element	corotTruss	1320	360	300	$A	$matTag
element	corotTruss	1321	361	301	$A	$matTag
element	corotTruss	1322	362	302	$A	$matTag
element	corotTruss	1323	363	303	$A	$matTag
element	corotTruss	1324	364	304	$A	$matTag
element	corotTruss	1325	365	305	$A	$matTag
element	corotTruss	1326	366	306	$A	$matTag
element	corotTruss	1327	367	307	$A	$matTag
element	corotTruss	1328	368	308	$A	$matTag
element	corotTruss	1329	369	309	$A	$matTag
element	corotTruss	1330	370	310	$A	$matTag
element	corotTruss	1331	371	311	$A	$matTag
element	corotTruss	1332	372	312	$A	$matTag
element	corotTruss	1333	373	313	$A	$matTag
element	corotTruss	1334	374	314	$A	$matTag
element	corotTruss	1335	375	315	$A	$matTag
element	corotTruss	1336	376	316	$A	$matTag
element	corotTruss	1337	377	317	$A	$matTag
element	corotTruss	1338	378	318	$A	$matTag
element	corotTruss	1339	379	319	$A	$matTag
element	corotTruss	1340	380	320	$A	$matTag
element	corotTruss	1341	381	321	$A	$matTag
element	corotTruss	1342	382	322	$A	$matTag
element	corotTruss	1343	383	323	$A	$matTag
element	corotTruss	1344	384	324	$A	$matTag
element	corotTruss	1345	385	325	$A	$matTag
element	corotTruss	1346	386	326	$A	$matTag
element	corotTruss	1347	387	327	$A	$matTag
element	corotTruss	1348	388	328	$A	$matTag
element	corotTruss	1349	389	329	$A	$matTag
element	corotTruss	1350	331	330	$A	$matTag
element	corotTruss	1351	332	331	$A	$matTag
element	corotTruss	1352	333	332	$A	$matTag
element	corotTruss	1353	334	333	$A	$matTag
element	corotTruss	1354	335	334	$A	$matTag
element	corotTruss	1355	336	335	$A	$matTag
element	corotTruss	1356	337	336	$A	$matTag
element	corotTruss	1357	338	337	$A	$matTag
element	corotTruss	1358	339	338	$A	$matTag
element	corotTruss	1359	340	339	$A	$matTag
element	corotTruss	1360	341	340	$A	$matTag
element	corotTruss	1361	342	341	$A	$matTag
element	corotTruss	1362	343	342	$A	$matTag
element	corotTruss	1363	344	343	$A	$matTag
element	corotTruss	1364	345	344	$A	$matTag
element	corotTruss	1365	346	345	$A	$matTag
element	corotTruss	1366	347	346	$A	$matTag
element	corotTruss	1367	348	347	$A	$matTag
element	corotTruss	1368	349	348	$A	$matTag
element	corotTruss	1369	350	349	$A	$matTag
element	corotTruss	1370	351	350	$A	$matTag
element	corotTruss	1371	352	351	$A	$matTag
element	corotTruss	1372	353	352	$A	$matTag
element	corotTruss	1373	354	353	$A	$matTag
element	corotTruss	1374	355	354	$A	$matTag
element	corotTruss	1375	356	355	$A	$matTag
element	corotTruss	1376	357	356	$A	$matTag
element	corotTruss	1377	358	357	$A	$matTag
element	corotTruss	1378	359	358	$A	$matTag
element	corotTruss	1379	330	359	$A	$matTag
element	corotTruss	1380	361	360	$A	$matTag
element	corotTruss	1381	362	361	$A	$matTag
element	corotTruss	1382	363	362	$A	$matTag
element	corotTruss	1383	364	363	$A	$matTag
element	corotTruss	1384	365	364	$A	$matTag
element	corotTruss	1385	366	365	$A	$matTag
element	corotTruss	1386	367	366	$A	$matTag
element	corotTruss	1387	368	367	$A	$matTag
element	corotTruss	1388	369	368	$A	$matTag
element	corotTruss	1389	370	369	$A	$matTag
element	corotTruss	1390	371	370	$A	$matTag
element	corotTruss	1391	372	371	$A	$matTag
element	corotTruss	1392	373	372	$A	$matTag
element	corotTruss	1393	374	373	$A	$matTag
element	corotTruss	1394	375	374	$A	$matTag
element	corotTruss	1395	376	375	$A	$matTag
element	corotTruss	1396	377	376	$A	$matTag
element	corotTruss	1397	378	377	$A	$matTag
element	corotTruss	1398	379	378	$A	$matTag
element	corotTruss	1399	380	379	$A	$matTag
element	corotTruss	1400	381	380	$A	$matTag
element	corotTruss	1401	382	381	$A	$matTag
element	corotTruss	1402	383	382	$A	$matTag
element	corotTruss	1403	384	383	$A	$matTag
element	corotTruss	1404	385	384	$A	$matTag
element	corotTruss	1405	386	385	$A	$matTag
element	corotTruss	1406	387	386	$A	$matTag
element	corotTruss	1407	388	387	$A	$matTag
element	corotTruss	1408	389	388	$A	$matTag
element	corotTruss	1409	360	389	$A	$matTag
element	corotTruss	1410	330	390	$A	$matTag
element	corotTruss	1411	331	390	$A	$matTag
element	corotTruss	1412	332	390	$A	$matTag
element	corotTruss	1413	333	390	$A	$matTag
element	corotTruss	1414	334	390	$A	$matTag
element	corotTruss	1415	335	390	$A	$matTag
element	corotTruss	1416	336	390	$A	$matTag
element	corotTruss	1417	337	390	$A	$matTag
element	corotTruss	1418	338	390	$A	$matTag
element	corotTruss	1419	339	390	$A	$matTag
element	corotTruss	1420	340	390	$A	$matTag
element	corotTruss	1421	341	390	$A	$matTag
element	corotTruss	1422	342	390	$A	$matTag
element	corotTruss	1423	343	390	$A	$matTag
element	corotTruss	1424	344	390	$A	$matTag
element	corotTruss	1425	345	390	$A	$matTag
element	corotTruss	1426	346	390	$A	$matTag
element	corotTruss	1427	347	390	$A	$matTag
element	corotTruss	1428	348	390	$A	$matTag
element	corotTruss	1429	349	390	$A	$matTag
element	corotTruss	1430	350	390	$A	$matTag
element	corotTruss	1431	351	390	$A	$matTag
element	corotTruss	1432	352	390	$A	$matTag
element	corotTruss	1433	353	390	$A	$matTag
element	corotTruss	1434	354	390	$A	$matTag
element	corotTruss	1435	355	390	$A	$matTag
element	corotTruss	1436	356	390	$A	$matTag
element	corotTruss	1437	357	390	$A	$matTag
element	corotTruss	1438	358	390	$A	$matTag
element	corotTruss	1439	359	390	$A	$matTag
element	corotTruss	1440	360	390	$A	$matTag
element	corotTruss	1441	361	390	$A	$matTag
element	corotTruss	1442	362	390	$A	$matTag
element	corotTruss	1443	363	390	$A	$matTag
element	corotTruss	1444	364	390	$A	$matTag
element	corotTruss	1445	365	390	$A	$matTag
element	corotTruss	1446	366	390	$A	$matTag
element	corotTruss	1447	367	390	$A	$matTag
element	corotTruss	1448	368	390	$A	$matTag
element	corotTruss	1449	369	390	$A	$matTag
element	corotTruss	1450	370	390	$A	$matTag
element	corotTruss	1451	371	390	$A	$matTag
element	corotTruss	1452	372	390	$A	$matTag
element	corotTruss	1453	373	390	$A	$matTag
element	corotTruss	1454	374	390	$A	$matTag
element	corotTruss	1455	375	390	$A	$matTag
element	corotTruss	1456	376	390	$A	$matTag
element	corotTruss	1457	377	390	$A	$matTag
element	corotTruss	1458	378	390	$A	$matTag
element	corotTruss	1459	379	390	$A	$matTag
element	corotTruss	1460	380	390	$A	$matTag
element	corotTruss	1461	381	390	$A	$matTag
element	corotTruss	1462	382	390	$A	$matTag
element	corotTruss	1463	383	390	$A	$matTag
element	corotTruss	1464	384	390	$A	$matTag
element	corotTruss	1465	385	390	$A	$matTag
element	corotTruss	1466	386	390	$A	$matTag
element	corotTruss	1467	387	390	$A	$matTag
element	corotTruss	1468	388	390	$A	$matTag
element	corotTruss	1469	389	390	$A	$matTag


 
 pattern Plain 1 Linear { 
load 390 0 0 -1e6
}


set tNode 390
recorder Node -file nodeDisp.txt  -time -node 390 -dof 1 2 3 disp

numberer Plain

system SparseSPD

constraints Transformation

test NormUnbalance 1e-5 8

algorithm Newton;

set arc_length 1.0;

#$type = 1 Minimum Residual Disp ;
#$type = 2 Normal Plain ;
#$type = 3 Update Normal Plain ;
#$type = 4 Cylindrical Arc-Length ;
set type 1;

integrator EQPath $arc_length $type



analysis Static;

puts "start"

set i 0;
set ret 0;
while {$i<165 && $ret==0} {
	incr i
	set ret [analyze 1];
}

puts "done"
