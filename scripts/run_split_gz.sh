#!/bin/bash



target_dir="/shares/grossniklaus.botinst.uzh/dkt/projects/meth1000/analysis/08_methcall_bedgraph/output/"

output_dir="/shares/grossniklaus.botinst.uzh/eharputluoglu/datasets_gz/"

error_log="$output_dir/error_log.txt"



cd "$target_dir" || exit 1



function process_files {

    # Count total number of files to process

    total_files=$(find . -maxdepth 2 -type f -name "*bismark.cov.gz" | wc -l)

    current_file=0



    # Loop through each subdirectory

    for dir in */; do

        cd "$dir" || continue



        mkdir -p "$output_dir/$dir"



        for file in *bismark.cov.gz; do

            base_name="${file%.bismark.cov.gz}"

            all_exist=true



            for chr in {1..5} Mt Pt; do

                if [ ! -f "$output_dir/$dir/${base_name}_chr${chr}.bismark.cov.gz" ]; then

                    all_exist=false

                    break

                fi

            done



            if [ "$all_exist" = false ]; then

                for chr in {1..5} Mt Pt; do

                    zcat "$file" | awk -v chr="$chr" -F '\t' '

                        $1 == chr && NF >= 6 {

                            if ($2 == $3) {  

                                $3 += 1        

                            }

                            printf "%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $5, $6

                        }' | gzip > "$output_dir/$dir/${base_name}_chr${chr}.bismark.cov.gz" 2>> "$error_log"

                done

            fi



            # Update progress percentage

            current_file=$((current_file + 1))

            percentage=$((current_file * 100 / total_files))

            printf "\rProgress: %3d%%" "$percentage"  # Print only the percentage

        done

        cd ..

    done

}





while true; do

    process_files

    

    # Check if all files have been processed

    if [ $current_file -eq $total_files ]; then

        break  # Exit the loop if all files are done

    else

        echo "Interrupted. Resuming..."

        sleep 5  # Wait for 5 seconds before restarting

    fi

done



echo ""

echo "Error log saved to: $error_log" 


