#!/bin/bash

target_dir="/shares/grossniklaus.botinst.uzh/dkt/projects/meth1000/analysis/08_methcall_bedgraph/output/"

output_dir="/shares/grossniklaus.botinst.uzh/eharputluoglu/datasets/"


cd "$target_dir" || exit 1


total_files=$(find . -maxdepth 2 -type f -name "*bismark.cov.gz" | wc -l)

current_file=0


# Loop through each subdirectory

for dir in */; do

    cd "$dir" || continue


    # Create sub dir

    mkdir -p "$output_dir/$dir"


    # Loop through each file matching the pattern *bismark.cov.gz

    for file in *bismark.cov.gz; do

        base_name="${file%.bismark.cov.gz}"

        all_exist=true


        # Check if all chromosome files already exist (with .gz extension)

        for chr in {1..5} Mt Pt; do 

            if [ ! -f "$output_dir/$dir/${base_name}_chr${chr}.bismark.cov.gz" ]; then

                all_exist=false

                break

            fi

        done


        # Process the file only if not all chromosome files exist

        if [ "$all_exist" = false ]; then

            for chr in {1..5} Mt Pt; do 

                zcat "$file" | awk -v chr="$chr" '

                    $1 == chr {

                        if ($2 == $3) {  

                            $3 += 1        

                        }

                        print

                    }' | gzip > "$output_dir/$dir/${base_name}_chr${chr}.bismark.cov.gz"  

            done

        else

            # Increment current_file even if skipping

            current_file=$((current_file + 1))

        fi


        # Update progress

        percentage=$((current_file * 100 / total_files))

        printf "\rProgress: [%-"$total_files"s] %3d%%" $(printf '=%*s' "$percentage" "") "$percentage"

    done

    cd ..

done

echo ""

