p=#plink dir
b=#beagle dir

# set variables 
mode=$1
data=$2
data_dir=$3
out_dir=$4
start_set=$5
end_set=$6
start_chr=$7
end_chr=$8
#exclude_list=$9
outfile=$9
ex_mode=$10
 
#set default values 
if [ "$data" == "" ]; then
     data=camh
fi
 
if [ "$data_dir" == "" ]; then
     data_dir=#target data
fi
 
if [ "$out_dir" == "" ]; then
    out_dir=#output dir
fi
 
if [ "$start_set" == "" ]; then
    start_set=1
fi
 
if [ "$end_set" == "" ]; then
    end_set=10
fi
 
if [ "$start_chr" == "" ]; then
    start_chr=6
fi
 
if [ "$end_chr" == "" ]; then
    end_chr=6
fi
 
#exclude="--exclude"
#if [ "$exclude_list" == "" ]; then
#    exclude=""
#fi
 
if [ "$outfile" == "" ]; then
    outfile=onekgimputation.out
fi
 
if [ "$ex_mode" == "" ]; then
    ex_mode="run"
fi
 
echo "1=$1 2=$2 3=$3 4=$4 5=$5 6=$6 7=$7 8=$8"
 
# ---------------------------------------------------------------------------------------
# Split by chromosome 
# ---------------------------------------------------------------------------------------
 
if [ "$mode" == "split.chr" ]; then
mkdir split_by_chr 
./split_bed_into_chromosomes.pl ${data_dir}/${data} split_by_chr/${data}.chr
fi 
 
# ---------------------------------------------------------------------------------------
# Generate a new dataset that includes only common SNPs in hapmap reference files 
# ---------------------------------------------------------------------------------------

if [ "$mode" == "split.chr.input" ]; then
    # number of bp in each split on chromosome & overlap of splits
	let split_size=3000000
    let overlap=750000
    let CHR=$start_chr CHR_NO=$end_chr
    
	#Loop through chromosomes 
    while [ ${CHR} -le ${CHR_NO} ]; do
        echo "CHR ${CHR} ... "
        
		#Make directory
        if [ ! -d chr${CHR} ]; then
            mkdir -p chr${CHR}
        fi
        
		# | <- pipes results from tail command to gawk function
		let start_pos=`head -1 region.${data}.chr${CHR}.bim | gawk '{ print \$4 }'`
        let end_pos=`tail -1 region.${data}.chr${CHR}.bim | gawk '{ print \$4 }'`
        let real_split_size=${split_size}-$overlap
        let chr_split_set=${end_pos}/${real_split_size}
        let is_last=(${end_pos}-${chr_split_set}*${split_size})
 
        if [ $is_last -ge $overlap ]; then
            let chr_split_set=${chr_split_set}+1                
        fi
 	   	# number of bp | number of chunks being split into
        echo "chr_size=${end_pos}  set_no=${chr_split_set}";       
		
        let s_i=1 s_e=${chr_split_set}
        while [ $s_i -le $s_e ]; do
            #* generate a list of splitted SNPs
            
			# if initial (s_i) equal to 1
			let start=(${start_pos}-1)
         echo ${start} 
			# if initial (s_i) not equal to 1 
			if [ ${s_i} -ne 1 ]; then
                start=${end}-${overlap}
            fi
            
			# calculate start and end bp positions 
            let start=${start}+1
            let end=${start}+${split_size}-1
             
            if [ ${s_i} -eq ${s_e} ]; then
                end=$end_pos
            fi
           
	    chmod 777 chr${CHR} 
		#write snps from bim file, between start and end positions, to new snps file
            gawk -v start=$start -v end=$end '$4>=start && $4<=end { print $2 }' region.${data}.chr${CHR}.bim > chr${CHR}/chr${CHR}.split.${s_i}.snps
            #number of snps within start-end range
			let snp_no=`wc -l chr${CHR}/chr${CHR}.split.${s_i}.snps | gawk '{ print \$1}'`
            echo "set ${s_i} : start_pos=${start} end_pos=${end}  snp_no=${snp_no}"
            
			
            if [ ${snp_no} -ge 0 ]; then
				#if less than 1000 snps, add next chunk to current one (+= split_size) and create new snps file
                if [ ${snp_no} -le 1000 ]; then
                    let end=${end}+${split_size}
                    gawk -v start=${start} -v end=${end} '$4>=start && $4<=end { print $2 }' region.${data}.chr${CHR}.bim > chr${CHR}/chr${CHR}.split.${s_i}.snps
                    let snp_no=`wc -l chr${CHR}/chr${CHR}.split.${s_i}.snps | gawk '{ print \$1}'`
                    echo "set ${s_i} (updated) : start_pos=${start} end_pos=${end}  snp_no=${snp_no}"
                fi
                
				# write "chr# | chunk # | starting bp # | ending bp #" to outfile
				echo "${CHR} ${s_i} ${start} ${end}" >> ${outfile}
                	
                #bsub -R "rusage[mem=3000]" -o ${out_dir}/chr${CHR}.log -e ${out_dir}/chr${CHR}.err -q hour \
                plink --bfile region.${data}.chr${CHR} --extract chr${CHR}/chr${CHR}.split.${s_i}.snps --recode-beagle --out chr${CHR}/split${s_i}
           fi
 
           let s_i=${s_i}+1
 
           if [ ${end} -ge ${end_pos} ]; then
            let s_i=${s_e}+100
           fi
 
        done
        let CHR=${CHR}+1
    done
fi
 
# ---------------------------------------------------------------------------------------
# 
# ---------------------------------------------------------------------------------------

if [ "$mode" == "gen.split.ref" ]; then
    list_file=$2
    let line=1
    let line_no=`wc -l $list_file | gawk '{ print \$1 }'`
     
    while [ $line -le $line_no ]; do
        # Extract from "split.chr.input" outfile
		let chr=`gawk -v line=$line 'NR==line { print \$1 }' $list_file`
        let set=`gawk -v line=$line 'NR==line { print \$2 }' $list_file`
        let start=`gawk -v line=$line 'NR==line { print \$3 }' $list_file`
        let end=`gawk -v line=$line 'NR==line { print \$4 }' $list_file`
 
        echo "chr=${chr} set=${set} start=${start} end=${end}"
		
		#create reference panel with no header
 	   	echo "$(tail -n +2 ref_panel_phased_header.bgl)" > ref_panel_phased.bgl
		
		#Extract snps belonging to each chunk, list ids and row number in .snps file; list all 4 fields in .markers file
        # extract reference SNPs based on chr pos
        # hap file: marker pos a1 a2
        gawk -v start=$start -v end=$end '$2>=start && $2<=end { print $1,NR }' markers_file.markers > chr${chr}/ref.chr${chr}.split.${set}.snps
        gawk -v start=$start -v end=$end '$2>=start && $2<=end { print $0 }' markers_file.markers > chr${chr}/ref.chr${chr}.split.${set}.begl.markers
        
		# starting and ending snp position
		let start_snp_line=`gawk 'NR==1 { print \$2}' chr${chr}/ref.chr$chr.split.$set.snps`
        let end_snp_line=`tail -1 chr${chr}/ref.chr$chr.split.$set.snps | gawk ' { print \$2}' `
 
        echo "start_snp_list=${start_snp_line} end_snp_line=${end_snp_line}"
 
        # extract phase information for reference SNPs 
        # phased file: M marker alleles
        gawk -v start=$start_snp_line -v end=$end_snp_line ' NR>=start && NR<=end { print $0 }' ref_panel_phased.bgl > chr${chr}/ref.chr${chr}.split.${set}.begl
        wc -l chr${chr}/ref.chr${chr}.split.${set}.begl
        wc -l chr${chr}/ref.chr${chr}.split.${set}.begl.markers
        wc -l chr${chr}/ref.chr${chr}.split.${set}.snps
 
        echo ""
 
        let line=${line}+1
    done
fi

# ---------------------------------------------------------------------------------------
# Impute 
# ---------------------------------------------------------------------------------------
if [ "$mode" == "step.impute" ]; then
    let CHR=${start_chr} CHR_NO=${end_chr}
	echo ${CHR}
    while [ $CHR -le $CHR_NO ]; do
			list_file=$2
	    	let split_start=1
	    	let split_end=`wc -l $list_file | gawk '{ print \$1 }'`
			
			#Make directory
			mkdir chr${CHR}_mod
			mkdir imputed			
			while [ $split_start -le $split_end ]; do
                    if [ -e "chr${CHR}/split${split_start}.chr-${CHR}.dat" ]; then
                        echo "CHR ${CHR} SPLIT ${split_start}: ${check_file}"
						
						#Extract C4 variant
						gawk '$2=="C4" { print $0 }' chr${CHR}/updated.ref.chr${CHR}.split.${split_start}.begl > chr${CHR}_mod/tmpC4.split.${split_start}.bgl
						#Remove C4 variant before strand flip
						gawk '$2!="C4" { print $0 }' chr${CHR}/updated.ref.chr${CHR}.split.${split_start}.begl > chr${CHR}_mod/trimmed.updated.ref.chr${CHR}.split.${split_start}.bgl
						#Markers file mod
						gawk '$1=="C4" { print $0 }' chr${CHR}/ref.chr${CHR}.split.${split_start}.begl.markers > chr${CHR}_mod/tmpC4.split.${split_start}.markers
						gawk '$1!="C4" { print $0 }' chr${CHR}/ref.chr${CHR}.split.${split_start}.begl.markers > chr${CHR}_mod/trimmed.updated.ref.chr${CHR}.split.${split_start}.markers						
						#Rename data files
						cp chr${CHR}/split${split_start}.chr-${CHR}.dat chr${CHR}_mod/split${split_start}.chr-${CHR}.bgl
						cp chr${CHR}/split${split_start}.chr-${CHR}.map chr${CHR}_mod/split${split_start}.chr-${CHR}.markers				
						#flip strands
						python check_strands/check_strands.py chr${CHR}_mod/trimmed.updated.ref.chr${CHR}.split.${split_start} chr${CHR}_mod/split${split_start}.chr-${CHR} chr${CHR}_mod/mod.split${split_start}.chr-${CHR}
						
						#remove blank files
						size=`ls -l chr${CHR}_mod/tmpC4.split.${split_start}.bgl|awk '{print $5}'`
						if [ $size -eq 0 ]
						  then
						    rm chr${CHR}_mod/tmpC4.split.${split_start}.bgl
						  fi
						size=`ls -l chr${CHR}_mod/tmpC4.split.${split_start}.markers|awk '{print $5}'`
						if [ $size -eq 0 ]
						  then
						    rm chr${CHR}_mod/tmpC4.split.${split_start}.markers
						  fi
						
						#put back C4
						if [ -e "chr${CHR}_mod/tmpC4.split.${split_start}.bgl" ]; then
							echo "Inserting C4" 
							gawk '$2>31949834 { print NR,$0 }' chr${CHR}_mod/mod.split${split_start}.chr-${CHR}.markers > chr${CHR}_mod/tmp_num.txt
							let NR=`gawk 'NR==1 { print $1 }' chr${CHR}_mod/tmp_num.txt` 
							line=`cat chr${CHR}_mod/tmpC4.split.${split_start}.markers`
							ex -sc "${NR}i|${line}" -cx chr${CHR}_mod/mod.split${split_start}.chr-${CHR}.markers
							join -1 2 -2 2 -o 1.2 1.3 <(sort -k2 chr${CHR}_mod/tmp_num.txt) <(sort -k2 chr${CHR}_mod/trimmed.updated.ref.chr${CHR}.split.${split_start}_mod.bgl) > chr${CHR}_mod/tmp_joined.txt
							sort -k2 chr${CHR}_mod/tmp_joined.txt > chr${CHR}_mod/tmp_joined_ordered.txt
							SNP=`gawk 'NR==1 { print $1 }' chr${CHR}_mod/tmp_joined_ordered.txt`
							let NR=`gawk -v snp=${SNP} '$2==snp { print NR }' chr${CHR}_mod/trimmed.updated.ref.chr${CHR}.split.${split_start}_mod.bgl`
							line=`cat chr${CHR}_mod/tmpC4.split.${split_start}.bgl`
							ex -sc "${NR}i|${line}" -cx chr${CHR}_mod/trimmed.updated.ref.chr${CHR}.split.${split_start}_mod.bgl
						fi
							#bsub -R "rusage[mem=6000]" -o ${out_dir}/chr${CHR}/begl.out -e ${out_dir}/chr${CHR}/begl.err \
                        
						#move files
						mv chr${CHR}_mod/split${split_start}.chr-${CHR}_mod.bgl imputed/split${split_start}.chr-${CHR}_mod.bgl
						mv chr${CHR}_mod/trimmed.updated.ref.chr${CHR}.split.${split_start}_mod.bgl imputed/trimmed.updated.ref.chr${CHR}.split.${split_start}_mod.bgl
						
						#Either use combined markers file (markers from both target data and reference panel)
						#mv chr${CHR}_mod/mod.split${split_start}.chr-${CHR}.markers imputed/mod.split${split_start}.chr-${CHR}.markers
						
						#or use only markers file corresponding to reference panel
						mv chr${CHR}/ref.chr${CHR}.split.${split_start}.begl.markers imputed/ref.chr${CHR}.split.${split_start}.begl.markers
						
						#Impute
						java -Xmx1000m -jar -Djava.io.tmpdir=imputed \
                        ${b} \
                        unphased=imputed/split${split_start}.chr-${CHR}_mod.bgl \
                        phased=imputed/trimmed.updated.ref.chr${CHR}.split.${split_start}_mod.bgl \
                        markers=imputed/ref.chr${CHR}.split.${split_start}.begl.markers  \
                        missing=0 out=imputed/bgl.chr${CHR}.split${split_start}
                    fi
            let split_start=${split_start}+1
            done
    let CHR=${CHR}+1
    done
fi 
