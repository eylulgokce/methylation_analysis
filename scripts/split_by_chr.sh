#!/bin/bash


target_dir="/shares/grossniklaus.botinst.uzh/dkt/projects/meth1000/analysis/08_methcall_bedgraph/output/"

output_dir="/shares/grossniklaus.botinst.uzh/eharputluoglu/datasets/"



cd "$target_dir" || exit 1


total_files=$(find . -maxdepth 2 -type f -name "*bismark.cov.gz" | wc -l)

current_file=0

for dir in */; do

    cd "$dir" || continue

    mkdir -p "$output_dir/$dir"


    for file in *bismark.cov.gz; do

        base_name="${file%.bismark.cov.gz}"

        all_exist=true  # if all output files exist


        # Check if all chromosome files already exist

        for chr in {1..5}; do

            if [ ! -f "$output_dir/$dir/${base_name}_chr${chr}.bismark.cov" ]; then

                all_exist=false

                break

            fi

        done


        # Process the file only if not all chromosome files exist
        if [ "$all_exist" = false ]; then
            for chr in {1..5}; do
                zcat "$file" | awk -v chr="$chr" '
                    $1 == chr {

                        if ($2 == $3) {  

                            $3 += 1        
                        }
                        print
                    }' | gzip > "$output_dir/$dir/${base_name}_chr${chr}.bismark.cov"
            done
        else
            current_file=$((current_file + 1))

        fi

        percentage=$((current_file * 100 / total_files))
        printf "\rProgress: [%-"$total_files"s] %3d%%" $(printf '=%*s' "$percentage" "") "$percentage"
    done

    cd ..

done

echo "" 

