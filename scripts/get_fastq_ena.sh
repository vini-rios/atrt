#!/bin/bash

usage(){
    echo "Usage: $0 <file_with_accession> [OPTIONS]"
    echo ""
    echo "Download files in parallel using wget from a file with Accession codes or a single Accession code."
    echo "Will download all FASTQ files associated with the provided Accessions."
    echo "Checks file integrity using MD5 checksums"
    echo ""
    echo "Requires wget, curl, xargs, md5sum"
    echo ""
    echo "Options:"
    echo ""
    echo "-h, --help              Show help on how to use this script"
    echo "-n, -p, --parallel      Number of parallel downloads (default: 10)"
    echo "-o, --outputdir         Path of the output directory (default: current directory)"
    echo "-r, --retry             Automatically retry failed downloads, up to 5 in parallel (default: no)"
    echo "-l, --list              Output a list of all files to be downloaded and their MD5s, then exit (default: no)"
    echo "                        The list won't be downloaded, just listed. Can be used as input file"
    echo ""
    echo "Input must be an single Accession code or a file with Accessions codes, one per line:"
    echo "SRR6413149"
    echo "SRR6413150"
    echo "SRR6413151"
    echo "..."
    echo ""
    echo "Examples:"
    echo "  $0 accession_list.txt -n 20 -o ./downloads/"
    echo "  $0 -n 6 SRR6413149 --outputdir /path/to/data/  2>log.txt  # Redirect errors to log.txt"
    echo "  $0 -o ./data/raw accession_list.txt  # Input can be anywhere"
    echo "  $0 PRJEB3251  # Uses defaults: 10 parallel downloads, current directory"
    exit 1
}

# Default values
input_file=""
parallel_jobs=10
output_dir="."
parallel_option_used=false
output_option_used=false
retry_option_used=false
list_only=false

# script with no arguments returns help
if [[ $# -eq 0 ]]; then
    usage
    exit 1
fi

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help) usage
        ;;
        -n|-p|--parallel)
            if [[ -n "$2" ]] && [[ ${2:0:1} != "-" ]]; then #If next argument is not empty and not another option (1st char -)
                if [[ "$parallel_option_used" = true ]]; then
                    echo "Warning: Parallel option specified multiple times, using latest value: $2" >&2
                fi
                if [[ "$2" =~ ^[1-9][0-9]*$ ]]; then # Check for positive integer
                    parallel_jobs="$2"
                    parallel_option_used=true
                    shift 2
                else
                    echo "Error: --parallel must be a positive integer" >&2
                    exit 1
                fi
            else
                echo "Error: Option $1 requires a positive integer, omit it for default 10" >&2
                exit 1
            fi
        ;;
        -o|--outputdir)
            if [[ -n "$2" ]] && [[ ${2:0:1} != "-" ]]; then
                if [[ "$output_option_used" = true ]]; then
                    echo "Warning: Output directory option specified multiple times, using latest value: $2" >&2
                fi
                output_option_used=true
                output_dir="$2"
                if [[ "${output_dir}" != */ ]]; then #add / if not present
                    output_dir="${output_dir}/" 
                fi
                shift 2
            else
                echo "Error: Option $1 requires an path, omit for current " >&2
                exit 1
            fi
        ;;
        -r|--retry)
            if [[ "$retry_option_used" = true ]]; then
                echo "Warning: Retry option specified multiple times, ignoring duplicates." >&2
            fi
            retry_option_used=true
            shift
        ;;
        -l|--list)
            echo "Listing all files to be downloaded and their MD5s, then exiting."
            list_only=true
            shift
        ;;
        -*) # Unknown option starting with -
            echo "Error: Unknown option $1" >&2
            echo "Use -h or --help for usage information."
            exit 1
        ;;
        *) # Positional argument (input file or accession code)
            if [[ -z "$input_file" ]]; then #If input_file is empty
                input_file="$1"
                shift
            else
                echo "Error: Multiple inputs provided. Only one accession code or file with multiples is allowed." >&2
                echo "For output directory use the -o or --outputdir option." >&2
                echo "First: $input_file"
                echo "Second: $1"
                exit 1
            fi
            ;;
    esac
done

# Validate that input_file was provided
if [[ -z "$input_file" ]]; then
    echo "Error: Input file or Accession code is required." >&2
    echo "Use -h or --help for usage information."
    exit 1
fi

# Check if input is a single accession code or a file
original_input="$input_file"
if [[ -f "$input_file" ]]; then
    echo "Input is a file: $input_file"
else
    echo "Input is a single Accession code: $input_file"
    temp_file=$(mktemp)
    echo "$input_file" > "$temp_file"
    input_file="$temp_file"
fi  

# Create output directory if it doesn't exist
if [[ ! -d "$output_dir" ]]; then
    echo "Creating output directory: $output_dir"
    mkdir -p "$output_dir" || {
        echo "Error: Failed to create output directory '$output_dir'" >&2
        exit 1
    }
fi

# Temporary files
ena_list=$(mktemp)
final_list=$(mktemp)
failed_list=$(mktemp)

# Read each line from input file and fetch associated FASTQ URLs and MD5s from ENA
while read -r line; do
    wget --quiet --output-document="$ena_list" "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$line&result=read_run&fields=run_accession,fastq_md5,fastq_ftp,library_layout"
    tail -n +2 $ena_list | while IFS=$'\t' read -r run_accession fastq_md5 fastq_ftp library_layout; do
        if [[ "$library_layout" == "PAIRED" ]]; then
            ftp_1="${fastq_ftp%%;*}" # %% removes longest match from end, everything after ;
            ftp_2="${fastq_ftp##*;}" # ## removes longest match from start, everything before ;
            md5_1="${fastq_md5%%;*}"
            md5_2="${fastq_md5##*;}"
            
            echo -e "$run_accession\t$ftp_1\t$md5_1\t$line" >> "$final_list"
            echo -e "$run_accession\t$ftp_2\t$md5_2\t$line" >> "$final_list"
        else
            md5_1="$fastq_md5"
            ftp_1="$fastq_ftp"
            echo -e "$run_accession\t$ftp_1\t$md5_1\t$line" >> "$final_list"
        fi
    done
done < "$input_file"

if [[ "$list_only" == true ]]; then
    cat $final_list > $output_dir/$original_input.files_to_download.tsv
    echo "List of files to be downloaded and their MD5s saved to: $output_dir$original_input.files_to_download.tsv"
    echo "Format: run_accession<TAB>ftp_link<TAB>md5<TAB>input_accession"
    exit 0
fi


echo ""
echo "Total files to download: $(wc -l < "$final_list")"
echo ""

download_fastq(){
    local run_accession="$1"
    local url="$2"
    local md5="$3"
    local original_accession="$4"
    local filename=$(basename "$url")
    echo "Downloading: $filename"
    wget --timestamping --continue --quiet --timeout=30 --tries=3 -O "${output_dir}/$4__${filename}" "ftp://$url"
    if [[ $? -eq 0 ]]; then
        echo "Download completed: $filename"
        local downloaded_file="${output_dir}/$4__${filename}"
        local calculated_md5=$(md5sum "$downloaded_file" | cut -d' ' -f1)
        if [[ "$calculated_md5" == "$md5" ]]; then
            echo "✓ MD5 verified: $filename"
        else
            echo "✗ MD5 mismatch for $filename (expected: $md5, got: $calculated_md5)" >&2
            echo -e "$run_accession\t$url\t$md5\t$original_accession" >> "$failed_list"
            rm -f "$downloaded_file"
            echo "Deleted corrupted file: $filename"
        fi
    else
        echo "Failed to download: $filename" >&2
        echo -e "$run_accession\t$url\t$md5\t$original_accession" >> "$failed_list"
    fi
}
export -f download_fastq # Export function for use in subshells
export output_dir # Export output_dir for use in subshells

# Start parallel downloads
cat "$final_list" | awk 'BEGIN{OFS=" "}{print $1,$2,$3,$4}' | xargs -P "$parallel_jobs" -n 4 bash -c 'download_fastq "$@"' _ # n 4 means 4 arguments per command, _ is a placeholder for $0


# Retry failed downloads if the option was used
failed_again=$(mktemp)
if [[ "$retry_option_used" == true ]]; then
    if [[ -s "$failed_urls" ]]; then # Check if file is not empty
        echo "Retrying failed downloads..."
        echo "Total failed downloads to retry: $(wc -l < "$failed_urls")"
        cat "$failed_urls" | awk 'BEGIN{OFS=" "}{print $1,$2,$3,$4}' | xargs -P 5 -n 4 bash -c 'download_fastq "$@"' _ # Retry with up to 5 in parallel
        if [[ -s "$failed_again" ]]; then
            echo "Some downloads still failed after retrying:" >&2
            cat "$failed_again" >&2
            exit 1
        else
            echo "All failed downloads succeeded on retry."
        fi
    else
        echo "No failed downloads to retry."
    fi
else
    if [[ -s "$failed_urls" ]]; then
        echo "Some downloads failed:" >&2
        cat "$failed_urls" >&2
        exit 1
    else
        echo "All downloads completed successfully."
    fi
fi

# Clean up temporary files
trap 'rm -f "$ena_list" "$final_list" "$failed_list" "failed_again" 2>/dev/null' EXIT