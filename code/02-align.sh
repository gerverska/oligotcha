# Define default parameters ####
threads=$(nproc --all)
evalue=1000
pid=50
word=10
targets=10000

# Define help message ####
usage=$(echo "-t number of threads to use
-e minimum expect value for hit
-p minimum percent identity for hit
-w minimum word size for hit
-m maximum number of hits per query
-h show help text")

# Update parameters if user supplied different arguments ####
while getopts ":t:e:p:w:m:h" option; do
    case $option in
        t) threads="$OPTARG";;
        e) evalue="$OPTARG";;
        p) pid="$OPTARG";;
        w) word="$OPTARG";;
        m) targets="$OPTARG";;
        h)
            echo "$usage"
            exit
            ;;
        \?)
            echo "Invalid option -$OPTARG" >&2
            echo "$usage" >&2
            exit 1
        ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            echo "$usage" >&2
            exit 1
        ;;
    esac
    
    case $OPTARG in
        -*)
            echo "Option $opt needs a valid argument"
            echo "$usage" >&2
            exit 1
        ;;
    esac
done

# Set input, output, log, and temp directories ####
in='01-design'
out='02-align'
db="${out}/db"
logs="${out}/logs"
rm -r $out
mkdir -p $db $logs
touch $out/README.md

# Decompress the reference FASTA ####
gzip -kd data/*silva-lsu-fungi.fa.gz

# Make the fungal database ####
makeblastdb -in data/*silva-lsu-fungi.fa -dbtype nucl \
-title "SILVA LSU Fungi" -out $db/silva-lsu-fungi \
-logfile $logs/silva-lsu-fungi.txt

# Align the blocking oligo sequence to the same database ####
blastn -query $in/blockers.fa -db $db/silva-lsu-fungi \
-dust no -evalue $evalue -perc_identity $pid -word_size $word -max_target_seqs $targets \
-num_threads $threads -outfmt '6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
-out $out/blockers.txt

# Remove the decompressed FASTA ####
rm data/*silva-lsu-fungi.fa
