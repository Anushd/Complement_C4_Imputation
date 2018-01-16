#set working directory
cd ~/anush/beagle_dir

for data in camh cidar-va gap imh2 mcic mgh-phrs-landr-nefs pages-con1 tcd-nuig1 umcu1 umcu3 catie cogs-uk imh1 imh3 mfs pages pages-con2 tcd-nuig2 umcu2
do
	#remove pre-existing output file
	if [ -e "onekgimputation.out" ]; then
	rm onekgimputation.out
	fi

	./imputation.sh 'split.chr' "${data}"
	./extract_from_bim_2.sh "${data}"
	./imputation.sh 'split.chr.input' "${data}"
	./imputation.sh 'gen.split.ref' 'onekgimputation.out'

	#add header to bgl files
	START=1
	file=onekgimputation.out
	END=`wc -l $file | gawk '{ print \$1 }'`
	while [[ $START -le $END ]]
	do
		awk 1 header.bgl chr6/ref.chr6.split.${START}.begl > chr6/updated.ref.chr6.split.${START}.begl
		((START = START + 1))
	done

	./imputation.sh 'step.impute' 'onekgimputation.out'

	#clean up
	rm -r chr6
	rm -r chr6_mod
	rm -r split_by_chr
	rm ref_panel_phased.bgl
	rm region.${data}*
	rm onekgimputation.out
	
	#rename imputed folder
	mv imputed imputed_${data}
done