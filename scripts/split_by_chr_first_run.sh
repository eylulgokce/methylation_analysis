#!/bin/bash


target_dir="/shares/grossniklaus.botinst.uzh/dkt/projects/meth1000/analysis/08_methcall_bedgraph/output/"
output_dir="/shares/grossniklaus.botinst.uzh/eharputluoglu/datasets/"

cd "$target_dir" || exit 1  

# Count total number of files to for progress bar
total_files=$(find . -maxdepth 2 -type f -name "*bismark.cov.gz" | wc -l)
current_file=0

# Loop through each subdirectory

for dir in */; do

    # Enter the subdirectory

    cd "$dir" || continue      


    # Create sub dir
    mkdir -p "$output_dir/$dir"

    # Loop through each file matching the pattern *bismark.cov.gz

    for file in *bismark.cov.gz; do

        # Extrac the base filenam

        base_name="${file%.bismark.cov.gz}"

        for chr in {1..5}; do
         
            zcat "$file" | awk -v chr="$chr" ' 
                $1 == chr {
                    if ($2 == $3) {   # check if start and end are same
                        $3 += 1        # increase end by 1
                    }
                    print
                }
            ' | gzip > "$output_dir/$dir/${base_name}_chr${chr}.bismark.cov"
        done



	# Update progress
        current_file=$((current_file + 1))
        percentage=$((current_file * 100 / total_files))
        printf "\rProgress: [%-"$total_files"s] %3d%%" $(printf '=%*s' "$percentage" "") "$percentage"

    done



    # Return to the parent directory

    cd ..

done

echo ""
