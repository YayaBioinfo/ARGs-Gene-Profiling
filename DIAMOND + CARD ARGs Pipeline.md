It performs **paired-end ARGs detection using DIAMOND against the CARD database**. It sets key thresholds for high-confidence identification: **minimum identity (`MIN_IDENTITY`) of 80%**, **minimum alignment length (`MIN_AMINO_ACIDS`) of 25 amino acids** (equivalent to 75 nucleotides, following ARGs-OAP recommendations), and an **E-value threshold (`EVALUE`) of 1e-7**. The pipeline automatically initializes directories, checks dependencies, validates the DIAMOND database, processes each sample’s paired-end files, runs DIAMOND BLASTX for both mates, combines results, validates outputs, and generates a detailed summary with total hits, high-confidence hits, average alignment length, and average identity.

```bash
#!/bin/bash

# ==================== CONFIGURATION ====================
readonly INPUT_DIR="/data/work/alignment_results/genome"
readonly DIAMOND_DB="/data/work/diamond_db/card_diamond"
readonly OUTPUT_DIR="/data/work/diamond_results"
readonly THREADS=8
readonly MIN_IDENTITY=80
readonly MIN_AMINO_ACIDS=25    # ARGs-OAP requirement: 25 aa = 75 nt
readonly EVALUE="1e-7"
readonly LOG_FILE="${OUTPUT_DIR}/analysis_log_$(date +%Y%m%d_%H%M%S).txt"

# ==================== INITIALIZATION ====================
init_analysis() {
    echo "=== 🚀 DIAMOND + CARD PAIRED-END ANALYSIS STARTED ==="
    echo "📅 Timestamp: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "📂 Input directory: $INPUT_DIR"
    echo "📁 Output directory: $OUTPUT_DIR"
    echo "🧬 DIAMOND DB: ${DIAMOND_DB}.dmnd"
    echo "🧵 Threads: $THREADS"
    echo "🎯 Min identity: ${MIN_IDENTITY}%"
    echo "📏 Min alignment: ${MIN_AMINO_ACIDS} aa (75 nt)"
    echo "🔬 E-value threshold: $EVALUE"
    echo "🔗 Mode: Paired-end (mate1 + mate2)"
    echo "=========================================="

    mkdir -p "$OUTPUT_DIR"

    # Logging
    exec > >(tee -a "$LOG_FILE") 2>&1
}

# ==================== DEPENDENCY CHECK ====================
check_dependencies() {
    local dependencies=("diamond")
    local missing=()

    for dep in "${dependencies[@]}"; do
        if ! command -v "$dep" &> /dev/null; then
            missing+=("$dep")
        fi
    done

    if [[ ${#missing[@]} -gt 0 ]]; then
        echo "❌ Missing dependencies: ${missing[*]}"
        echo "💡 Installing via mamba..."
        mamba install -c bioconda "${missing[@]}" -y
    else
        echo "✅ All dependencies satisfied"
    fi
}

# ==================== DATABASE VALIDATION ====================
validate_diamond_database() {
    if [[ -f "${DIAMOND_DB}.dmnd" ]]; then
        echo "✅ DIAMOND database found: ${DIAMOND_DB}.dmnd"
        return 0
    else
        echo "❌ ERROR: DIAMOND database missing: ${DIAMOND_DB}.dmnd"
        return 1
    fi
}

# ==================== SAMPLE PROCESSING ====================
process_sample() {
    local mate1_file="$1"
    local sample_base
    sample_base=$(basename "$mate1_file" _genome_Unmapped.out.mate1)
    local mate2_file="${INPUT_DIR}/${sample_base}_genome_Unmapped.out.mate2"

    local output_file="${OUTPUT_DIR}/${sample_base}_diamond_card.tsv"
    local temp_dir="${OUTPUT_DIR}/temp_${sample_base}"

    echo ""
    echo "🧬 Processing sample: $sample_base"
    echo "📄 Mate 1: $(basename "$mate1_file")"
    echo "📄 Mate 2: $(basename "$mate2_file")"

    if [[ ! -f "$mate1_file" || ! -f "$mate2_file" ]]; then
        echo "❌ ERROR: Missing paired-end files for $sample_base"
        return 1
    fi

    mkdir -p "$temp_dir"

    local start_time=$(date +%s)

    # Analyze mates
    analyze_mate "$mate1_file" "${temp_dir}/mate1.tsv" "mate1" || return 1
    analyze_mate "$mate2_file" "${temp_dir}/mate2.tsv" "mate2" || return 1

    # Combine results
    cat "${temp_dir}/mate1.tsv" "${temp_dir}/mate2.tsv" > "$output_file"

    rm -rf "$temp_dir"

    local end_time=$(date +%s)
    local elapsed=$((end_time - start_time))

    validate_results "$output_file" "$sample_base" "$elapsed"
}

# ==================== MATE ANALYSIS ====================
analyze_mate() {
    local input_file="$1"
    local output_file="$2"
    local mate="$3"

    echo "   🔬 Running DIAMOND on $mate..."

    diamond blastx \
        --query "$input_file" \
        --db "$DIAMOND_DB" \
        --out "$output_file" \
        --outfmt 6 \
        --threads "$THREADS" \
        --evalue "$EVALUE" \
        --id "$MIN_IDENTITY" \
        --query-cover "$MIN_AMINO_ACIDS" \
        --max-target-seqs 1 \
        --quiet

    return $?
}

# ==================== RESULT VALIDATION ====================
validate_results() {
    local file="$1"
    local sample="$2"
    local time="$3"

    if [[ ! -f "$file" ]]; then
        echo "❌ ERROR: Output file not created for $sample"
        return 1
    fi

    local total_hits
    total_hits=$(wc -l < "$file")

    local high_conf
    high_conf=$(awk -v id="$MIN_IDENTITY" '$3 >= id {c++} END {print c+0}' "$file")

    local avg_len
    avg_len=$(awk -F'\t' '{sum+=$4} END {print (NR?sum/NR:0)}' "$file")

    local avg_id
    avg_id=$(awk -F'\t' '{sum+=$3} END {print (NR?sum/NR:0)}' "$file")

    echo "✅ SUCCESS: $sample"
    echo "   ⏱️ Time: ${time}s"
    echo "   📊 Total hits: $total_hits"
    echo "   🎯 High-confidence hits: $high_conf"
    echo "   📏 Avg alignment length: $avg_len"
    echo "   🎯 Avg identity: $avg_id"

    return 0
}

# ==================== SUMMARY ====================
generate_summary() {
    local total="$1"
    local success="$2"

    echo ""
    echo "🎉 ANALYSIS COMPLETED"
    echo "=========================================="
    echo "📊 SUMMARY REPORT"
    echo "   Total samples: $total"
    echo "   Successful: $success"
    echo "   Failed: $((total - success))"
    echo "   🔥 Success rate: $((success * 100 / total))%"
    echo ""
    echo "📁 Results directory: $OUTPUT_DIR"
    echo "📋 Log file: $LOG_FILE"
    echo "=========================================="

    create_detailed_summary
}

create_detailed_summary() {
    local summary="${OUTPUT_DIR}/analysis_summary_$(date +%Y%m%d).csv"
    echo "Sample,Total_Hits,High_Confidence_Hits,Avg_Length,Avg_Identity" > "$summary"

    for file in ${OUTPUT_DIR}/*_diamond_card.tsv; do
        [[ -f "$file" ]] || continue

        local sample
        sample=$(basename "$file" _diamond_card.tsv)

        local hits
        hits=$(wc -l < "$file")

        local high
        high=$(awk -v id="$MIN_IDENTITY" '$3 >= id {c++} END {print c+0}' "$file")

        local len
        len=$(awk -F'\t' '{sum+=$4} END {print (NR?sum/NR:0)}' "$file")

        local ident
        ident=$(awk -F'\t' '{sum+=$3} END {print (NR?sum/NR:0)}' "$file")

        echo "$sample,$hits,$high,$len,$ident" >> "$summary"
    done

    echo "📈 Detailed summary saved to: $summary"
}

# ==================== MAIN ====================
main() {
    init_analysis
    check_dependencies

    if ! validate_diamond_database; then
        echo "❌ Critical error: database validation failed."
        exit 1
    fi

    local total=0
    local success=0

    echo ""
    echo "🚀 Starting paired-end analysis with ARGs-OAP parameters..."
    echo "----------------------------------------"

    for mate1 in ${INPUT_DIR}/*_genome_Unmapped.out.mate1; do
        ((total++))
        if process_sample "$mate1"; then
            ((success++))
        fi
        echo "---"
    done

    generate_summary "$total" "$success"
}

main "$@"
```
