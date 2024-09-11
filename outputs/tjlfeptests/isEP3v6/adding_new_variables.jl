using Printf
omegaGAM = [0.37630465555519, 0.375361034241556, 0.374380032336559, 0.37336088641571, 0.372302811806779, 0.371204973881369, 0.370066510887835, 0.368886503650228, 0.367664017421013, 0.366398038669372, 0.365087522842374, 0.36373136309873, 0.362328382666669, 0.360877357074407, 0.359890597878842, 0.359701190604852, 0.359496390666809, 0.359276256383357, 0.359040887898204, 0.358826076391722, 0.358688108632774, 0.358535385094891, 0.35836797206387, 0.358206963247139, 0.358081535209476, 0.357942093838046, 0.357788743760292, 0.357708083745107, 0.357622368112824, 0.357522355452742, 0.357466214762942, 0.357402707420741, 0.357326593676875, 0.357250573183121, 0.357160022215152, 0.357049311110812, 0.356920501379736, 0.356774019802091, 0.356595082932764, 0.356402785899601, 0.356220351558877, 0.356027318723833, 0.355832245715287, 0.355627570417849, 0.355441178366132, 0.35525496992423, 0.35508627058222, 0.35491608951371, 0.354746313530198, 0.354566051227472, 0.354442722365749, 0.354310490536241, 0.354175784055775, 0.354035535556911, 0.35389157792578, 0.353772686367281, 0.353654671284596, 0.353529569724462, 0.353395731885065, 0.353254077037527, 0.353149332954925, 0.353083141424206, 0.353093287119249, 0.353179858048903, 0.353337513242025, 0.353580961989492, 0.353913797737129, 0.354320420532363, 0.354804128917787, 0.355361467538128, 0.355986923871072, 0.356686738447008, 0.357463699726023, 0.358301144808805, 0.359196630437856, 0.360157289369343, 0.361173885786164, 0.362246300623102, 0.363377620479413, 0.364562189865451, 0.365796991590063, 0.367051564142604, 0.368303524829443, 0.369546607151744, 0.370780801243314, 0.372007537018915, 0.37322013079556, 0.374421274300664, 0.375611740032507, 0.376788344541803, 0.377955052499745, 0.379112083571969, 0.380259016613174, 0.381393151625539, 0.382519405094984, 0.383637716128822, 0.384746585880212, 0.385851198399325, 0.386947837002722, 0.38803680649957, 0.389120889044541, 0.390175365614123, 0.391184679555325, 0.392148953671739, 0.393066184005639, 0.393936857276485, 0.394759347376739, 0.395533585837782, 0.396260024066002, 0.39693769469877, 0.397570908485845, 0.398158507756068, 0.398696912103879, 0.399184313338574, 0.399622528130653, 0.400011211617096, 0.400350500509309, 0.400644572788768, 0.400885163799425, 0.40107172878629, 0.401202720684521, 0.401308839408675, 0.401421264282795, 0.40154231153316, 0.401674075582934, 0.401817871826159, 0.401959793475625, 0.402104988044211, 0.402257366027092, 0.402418740450103, 0.402592340140439, 0.402783103170248, 0.402992346558116, 0.403223588741344, 0.403474090819543, 0.40374697930223, 0.404049438781056, 0.404383760271533, 0.404754998954322, 0.405175641895724, 0.405634086130728, 0.406130956991703, 0.406662150955465, 0.407230127318191, 0.407840742516122, 0.408498181629081, 0.409196702730679, 0.409946398411876, 0.410745361097, 0.411587207860938, 0.412476045161064, 0.41340765401742, 0.414381960737918, 0.415401241365133, 0.416462802751122, 0.4175637458677, 0.41870533879622, 0.419885032609544, 0.421100358274037, 0.42235449065239, 0.423643162141398, 0.42495412855419, 0.42627523011523, 0.427604423635373, 0.428939751259235, 0.430278960820523, 0.431620579516189, 0.432969919503956, 0.43431823599528, 0.435665771929741, 0.437012076147528, 0.438355747155728, 0.439696838232402, 0.441038873145878, 0.442375838864776, 0.443706419067509, 0.445030391551356, 0.446361076520414, 0.447681929668527, 0.448988792039087, 0.450280131157404, 0.451263560105778, 0.452212813681002, 0.453123537985476, 0.453952995488533, 0.454701492262152, 0.455361020357383, 0.455850255033451, 0.456268587914944, 0.456602984464372, 0.456907577024117, 0.457235093488514, 0.45761250415146, 0.458118154695441, 0.458708823339313, 0.459421035121657, 0.460268975750747, 0.461168800185261, 0.462151550873789, 0.463182205866487, 0.464216522349101]
gammap = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
gammaE = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]



open("/home/towlej/.julia/dev/TJLF.jl/outputs/tjlfeptests/isEP3v6/input.EXPRO", "a") do file
    for (i, value) in enumerate(omegaGAM)
        # Write in the required format: EXPRO_omegaGAM i= value with specific spacing
        println(file, @sprintf("EXPRO_omegaGAM_%3d= %16.10f", i, value))
    end
    for (i, value) in enumerate(gammap)
        # Write in the required format: EXPRO_omegaGAM i= value with specific spacing
        println(file, @sprintf("EXPRO_gammap_%3d= %16.10f", i, value))
    end
    for (i, value) in enumerate(gammaE)
        # Write in the required format: EXPRO_omegaGAM i= value with specific spacing
        println(file, @sprintf("EXPRO_gammaE_%3d= %16.10f", i, value))
    end
end