FOR=ifort
DIR="./bin/"

all :
	$(FOR) dirgapf.f90 -o $(DIR)dirgapf.x
	$(FOR) fake_scf_kpath.f90 -o $(DIR)f_scf_kpath.x
	$(FOR) vasp_bands_kpath.f90 -o $(DIR)vbands_k.x
	$(FOR) vasp_bands_kpath_fake.f90 -o $(DIR)vbands_k_f.x
	$(FOR) vasp_bands_kpath_fake_procar.f90 -o $(DIR)vbands_k_f_p.x
	$(FOR) vasp_bands_kpath_fake_sp.f90 -o $(DIR)vbands_k_f_sp.x
	$(FOR) vasp_bands_kpath_scs.f90 -o $(DIR)vbands_k_scs.x
	$(FOR) vasp_dos_2-mod.f90 -o $(DIR)vaspdosj.x
	$(FOR) vasp_gap_vbm-finder.f90 -o $(DIR)gapf.x
	$(FOR) vasp_gap_vbm-finder-hse.f90 -o $(DIR)gapf_h.x
	$(FOR) vasp_gap_vbm-finder_sp.f90 -o $(DIR)gapf_sp.x
	$(FOR) vasp_gap_vbm-finder_sp-hse.f90 -o $(DIR)gapf_sp_h.x
	$(FOR) vasp_gap-finder-soc.f90 -o $(DIR)gaps.x
	$(FOR) vasp_gap-finder-soc-hse.f90 -o $(DIR)gaps_h.x
	$(FOR) vasp_gap_vbm-finder-sm.f90 -o $(DIR)gapf_sm.x
	$(FOR) rkmesh-POSCAR.f90 -o $(DIR)rkmesh.x
	$(FOR) ir-raman-spec.f90 -o $(DIR)ir-raman.x	







 
