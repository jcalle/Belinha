			#Jorge 2017 - Implementing error analysis to MHM Method
			add_subdirectory(MHM_Jorge)			
			#Jorge 2018 - Implementing error analysis from Paredes Thesis
			add_subdirectory(HDiv_Jorge)							
		#Testing generating ellipse from a vector of points by least square method
			add_subdirectory(EllipseGenerator)

			# Tutorial about hp-adaptive techniques
			add_subdirectory(Poisson3D_Fichera)
			add_subdirectory(Poisson3D_Shock)
			add_subdirectory(Poisson3D_AutoHP)

	add_subdirectory(Laplacian2D_GID)
	add_subdirectory(Elasticity2D_Wav_Gid)


	#Project developed by Olivier
	#add_subdirectory(OlivierRoussel)


#Projects to Conservation Laws
	add_subdirectory(ThesisPhD)
