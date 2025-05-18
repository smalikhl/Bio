
**Die übergeordneten Ziele der Hausaufgaben waren:**

1.  **Datenvorbereitung und Analyse für alle Proben:** Nicht nur eine Probe, sondern alle sechs (Noss.R1-R3, Pais.R1-R3).
2.  **Qualitätskontrolle:** Mit einem Werkzeug namens `seqstats` die Qualität der Daten vor und nach der Bearbeitung prüfen.
3.  **Ergebnisse verpacken:** Die wichtigen Ergebnisdateien für die Abgabe zusammenfassen.
4.  **Verständnis dokumentieren:** Die Schritte der Analyse in eigenen Worten erklären.

Dazu hast du drei Bash-Skripte verwendet. Ein Bash-Skript ist wie ein Kochrezept für den Computer, das ihm sagt, welche Befehle er in welcher Reihenfolge ausführen soll.

---

**Teil 1: Die PDF-Anleitungen – Unser Bauplan**

Bevor wir in die Skripte eintauchen, erinnern wir uns an die wichtigsten Vorgaben aus den PDF-Dokumenten (insbesondere "5 Implement an RNA-seq workflow" und die Folien "Practical class Introduction to Unix and RNA-seq analysis" sowie "Practical class RNA-seq data analysis").

*   **Unix-Grundlagen (PDF "Introduction to Unix...")**: Hier ging es darum, wie man sich im Terminal bewegt (`cd`), Dateien auflistet (`ls`), Verzeichnisse erstellt (`mkdir`), Dateien kopiert (`cp`), verschiebt/umbenennt (`mv`) und löscht (`rm`). Diese Befehle sind das Handwerkszeug. Auch das Komprimieren/Dekomprimieren (`gzip`/`gunzip`, `zip`/`unzip`) wurde erwähnt (siehe 1.8).
*   **RNA-Seq Workflow (PDF "5 Implement an RNA-seq workflow" & Folien "RNA-seq data analysis"):**
    *   **Ziel:** Eine RNA-Seq Pipeline durchführen (siehe Folie "Overview of a minimal RNA-seq pipeline" und "Implement an RNA-seq workflow").
    *   **Schritt 1: Adapter- und Qualitäts-Trimming (PDF 5.1, Folie "Adapter and quality trimming")**
        *   **Warum?** Sequenziergeräte fügen kleine technische DNA-Stücke hinzu (Adapter) und die Qualität der Sequenzierung kann an den Enden der Reads schlecht sein. Das muss weg!
        *   **Werkzeug:** `cutadapt` (oder `trim-galore`).
        *   **Wichtige Parameter:**
            *   Adaptersequenz: `"AGATCGGAAGAG"` (für beide Reads).
            *   Nur 3'-Adapter entfernen.
            *   Paired-End Verarbeitung (beide Reads eines Paares zusammen behandeln).
            *   Qualitäts-Cutoff: `20`.
            *   **Minimale Leselänge nach Trimmen: `80` Basen.** (Das war ein wichtiger Punkt!)
    *   **Schritt 2: Reads Mapping (PDF 5.2, Folie "Mapping")**
        *   **Warum?** Herausfinden, woher im Genom (unserer Referenzkarte) die getrimmten Reads ursprünglich stammen.
        *   **Werkzeug:** `Hisat2` (oder `STAR`).
        *   **Benötigt:** Einen Genom-Index (eine Art Inhaltsverzeichnis für das Genom, damit die Suche schneller geht). Für `Ahalleri.fa` war der Index `genome/Ahalleri`.
    *   **Schritt 3: Reads Counting (PDF 5.3, Folie "Counting")**
        *   **Warum?** Zählen, wie viele Reads auf jedes bekannte Gen gemappt wurden. Das gibt uns ein Maß für die Genaktivität.
        *   **Werkzeug:** `featureCounts` (oder `htseq-count`).
        *   **Benötigt:** Eine Annotationsdatei (`.gff` oder `.gtf`), die sagt, wo die Gene im Genom liegen (hier `genome/Ahalleri.gff`).
    *   **Conda (PDF 3.3, 3.4, Folie "Bioinformatics software ... CONDA"):** Ein Werkzeug, um andere Bioinformatik-Programme und ihre Abhängigkeiten sauber zu installieren und zu verwalten, oft in isolierten Umgebungen. Die PDF empfiehlt "One tool = one environment" oder zumindest isolierte Umgebungen. Du hast eine Gesamtumgebung `rna_seq_env` erstellt, was für die Ausführung der Skripte sehr praktisch ist.

---

**Teil 2: Deine Bash-Skripte und die Ergebnisse im Detail**

Du hast drei Skripte verwendet:
1.  `run_rna_seq_pipeline.sh` (für die Hauptanalyse)
2.  `run_seqstats_analysis.sh` (für die Qualitätskontrolle mit seqstats)
3.  `create_submission_zip.sh` (zum Verpacken der Ergebnisse)

**Skript 1: `run_rna_seq_pipeline.sh` – Die Hauptanalyse**

**Was dieses Skript tut (Hausaufgabe a):**
Es nimmt deine Rohdaten, bereitet sie auf und führt die Kernanalyse für alle 6 Probenpaare durch.

```bash
#!/bin/bash
# Dies ist die "Shebang"-Zeile. Sie sagt dem Betriebssystem: "Hallo, ich bin ein Bash-Skript!
# Bitte benutze das Programm /bin/bash, um mich auszuführen."

echo "=== STARTING RNA-SEQ PIPELINE (Homework Version) ==="
# echo gibt einfach den Text dahinter auf dem Bildschirm aus. Eine nette Startnachricht.
SECONDS=0
# SECONDS ist eine eingebaute Variable in Bash. Wenn man ihr 0 zuweist, fängt sie an,
# die Sekunden zu zählen. Am Ende können wir so sehen, wie lange das Skript gelaufen ist.

# --- 0. Verzeichnisse und Pfade definieren ---
# Hier legen wir fest, wo sich wichtige Dateien und Ordner befinden.
# Das macht das Skript übersichtlicher und leichter anzupassen.
BASE_DIR=$(pwd)
# pwd (print working directory) gibt das aktuelle Verzeichnis aus, in dem du das Skript gestartet hast.
# $(...) führt den Befehl darin aus und gibt das Ergebnis zurück.
# BASE_DIR speichert also den Pfad zu deinem Haupt-Arbeitsverzeichnis
# (z.B. /seminar/20250414_BioInf/doulos36/RNA_seq/hausaufgabe_rna_seq).

# --- KORREKTUR BASIEREND AUF DEINER STRUKTUR ---
# Hier wurde der Pfad zu den Rohdaten angepasst, basierend auf deiner Verzeichnisstruktur.
RAW_DIR_INPUT="${BASE_DIR}/raw_reads_all"
# RAW_DIR_INPUT ist das Verzeichnis, wo das Skript nach komprimierten Rohdaten sucht
# (z.B. *.fq.gz, *.fq.zip). In deinem Fall ist das jetzt dein "raw_reads_all"-Ordner.
RAW_DIR_PROCESSED="${BASE_DIR}/raw_reads_all"
# RAW_DIR_PROCESSED ist das Verzeichnis, in dem die dekomprimierten Rohdaten (.fq)
# landen (falls sie dekomprimiert werden müssen) UND von dem aus Cutadapt die .fq-Dateien liest.
# Da deine .fq-Dateien schon in raw_reads_all sind, ist es hier dasselbe Verzeichnis.

TRIM_DIR="${BASE_DIR}/trimmed_reads"       # Ordner für getrimmte Reads
MAP_DIR="${BASE_DIR}/mapped_reads"         # Ordner für Mapping-Ergebnisse (SAM/BAM)
COUNT_DIR="${BASE_DIR}/counts_all"         # Ordner für Read-Counts
GENOME_DIR="${BASE_DIR}/genome"            # Dein Ordner mit Genomdaten (Ahalleri.fa, Ahalleri.gff, Indexe)
GENOME_ANNOTATION="${GENOME_DIR}/Ahalleri.gff" # Die GFF-Datei, die sagt, wo Gene sind.
HISAT2_INDEX_PREFIX="${GENOME_DIR}/Ahalleri"   # Das Präfix für die Hisat2-Indexdateien.

mkdir -p "$RAW_DIR_PROCESSED" "$TRIM_DIR" "$MAP_DIR" "$COUNT_DIR"
# mkdir -p erstellt die genannten Ordner, falls sie noch nicht da sind.
# -p sorgt dafür, dass auch übergeordnete Ordner erstellt werden (hier nicht nötig)
# und dass es keine Fehlermeldung gibt, wenn die Ordner schon existieren.
# GENOME_DIR wird nicht erstellt, da es laut deiner ls-Ausgabe schon da ist.

# --- 1. Dekomprimieren der Rohdaten (Hausaufgabe a) ---
# Ziel: Sicherstellen, dass wir unkomprimierte .fq-Dateien für die Analyse haben.
echo "--- [1/5] Dekomprimiere Rohdaten ---"

echo "Dekprimiere .gz Dateien aus $RAW_DIR_INPUT..."
# Eine Schleife, die alle Dateien im Ordner RAW_DIR_INPUT durchgeht, die auf ".fq.gz" enden.
for gz_file in "$RAW_DIR_INPUT"/*.fq.gz; do
    # if [ -f "$gz_file" ]; then: Prüft, ob die gefundene .gz-Datei wirklich eine Datei ist.
    if [ -f "$gz_file" ]; then
        target_fq_name=$(basename "${gz_file%.gz}") # Extrahiert den Dateinamen ohne Pfad und ohne .gz
        target_fq_path="$RAW_DIR_PROCESSED/$target_fq_name" # Der Zielpfad für die .fq-Datei
        # if [ ! -f "$target_fq_path" ]; then: Prüft, ob die .fq-Datei noch NICHT existiert.
        if [ ! -f "$target_fq_path" ]; then
            echo "  Dekomprimiere $gz_file zu $target_fq_path"
            gunzip -c "$gz_file" > "$target_fq_path" # gunzip -c dekomprimiert und schreibt auf Bildschirm, > leitet es in Datei um.
        else
            # Deine Ausgabe zeigt diesen Fall:
            echo "  Zieldatei $target_fq_path existiert bereits, überspringe gunzip für $gz_file."
        fi
    else
        # Diese Warnung kam bei dir nicht, da deine .gz-Dateien in raw_reads_all waren.
        echo "  WARNUNG: Quelldatei für .gz-Muster ($gz_file) nicht in $RAW_DIR_INPUT gefunden oder ungültig."
    fi
done # Ende der .gz-Schleife

echo "Dekprimiere .zip Dateien aus $RAW_DIR_INPUT..."
# Analog für .zip-Dateien.
for zip_file in "$RAW_DIR_INPUT"/*.fq.zip; do
    if [ -f "$zip_file" ]; then
        target_fq_name=$(basename "${zip_file%.zip}")
        target_fq_path="$RAW_DIR_PROCESSED/$target_fq_name"
        if [ ! -f "$target_fq_path" ]; then
            echo "  Dekomprimiere $zip_file nach $RAW_DIR_PROCESSED"
            unzip -q -n "$zip_file" -d "$RAW_DIR_PROCESSED" # -q (quiet), -n (never overwrite), -d (directory)
        else
            # Deine Ausgabe zeigt diesen Fall:
            echo "  Zieldatei $target_fq_path (aus $zip_file) existiert bereits, unzip wird sie nicht überschreiben (wegen -n Option)."
        fi
    else
        echo "  WARNUNG: Quelldatei für .zip-Muster ($zip_file) nicht in $RAW_DIR_INPUT gefunden oder ungültig."
    fi
done
echo "--- Dekompression abgeschlossen ---"
# **Bedeutung für dich:** Das Skript hat korrekt erkannt, dass deine .fq-Dateien schon in raw_reads_all sind
# und hat sie nicht unnötig neu entpackt. Das ist gut!

# --- 2. RNA-Seq Analyse-Schleife für alle 6 Paare ---
# Jetzt kommt der Hauptteil: Die Analyse für jede Probe.
SAMPLES=("Noss.R1" "Noss.R2" "Noss.R3" "Pais.R1" "Pais.R2" "Pais.R3")
# Ein Array (eine Liste) mit den Namen deiner Proben-Präfixe.

# Wichtige Parameter für cutadapt, wie in der PDF gefordert:
ADAPTER_SEQ="AGATCGGAAGAG" # Die zu entfernende Adaptersequenz.
QUALITY_CUTOFF=20          # Basen mit Qualität unter 20 am Ende werden abgeschnitten.
MIN_LENGTH_AFTER_TRIM=80   # Reads, die nach dem Trimmen kürzer als 80 Basen sind, werden verworfen.

# Die for-Schleife geht jedes Element im SAMPLES-Array durch.
# In jedem Durchlauf ist SAMPLE_PREFIX z.B. "Noss.R1", dann "Noss.R2" usw.
for SAMPLE_PREFIX in "${SAMPLES[@]}"; do
    echo ""
    echo "-----------------------------------------------------"
    echo "--- Verarbeite Sample: $SAMPLE_PREFIX ---" # Zeigt an, welche Probe gerade dran ist.
    echo "-----------------------------------------------------"

    # Dateinamen für die aktuelle Probe zusammenbauen:
    R1_RAW_FQ="$RAW_DIR_PROCESSED/${SAMPLE_PREFIX}.1.sub.fq" # z.B. raw_reads_all/Noss.R1.1.sub.fq
    R2_RAW_FQ="$RAW_DIR_PROCESSED/${SAMPLE_PREFIX}.2.sub.fq" # z.B. raw_reads_all/Noss.R1.2.sub.fq

    R1_TRIMMED_FQ="$TRIM_DIR/${SAMPLE_PREFIX}.1.trimmed.fq"   # Zieldatei für getrimmte R1
    R2_TRIMMED_FQ="$TRIM_DIR/${SAMPLE_PREFIX}.2.trimmed.fq"   # Zieldatei für getrimmte R2
    CUTADAPT_LOG="$TRIM_DIR/${SAMPLE_PREFIX}.cutadapt.log"    # Logdatei für cutadapt

    SAM_OUT="$MAP_DIR/${SAMPLE_PREFIX}.aligned.sam"           # SAM-Datei vom Mapping
    SORTED_BAM_OUT="$MAP_DIR/${SAMPLE_PREFIX}.aligned.sorted.bam" # Sortierte BAM-Datei
    HISAT2_SUMMARY="$MAP_DIR/${SAMPLE_PREFIX}.hisat2.summary.txt" # Hisat2 Statistik
    
    COUNT_FILE_PREFIX="$COUNT_DIR/${SAMPLE_PREFIX}" # Präfix für featureCounts-Dateien

    # Sicherheitscheck: Sind die .fq-Dateien da?
    if [ ! -f "$R1_RAW_FQ" ] || [ ! -f "$R2_RAW_FQ" ]; then
        echo "  FEHLER: Rohdatendateien (.fq) für $SAMPLE_PREFIX nicht gefunden in $RAW_DIR_PROCESSED ($R1_RAW_FQ oder $R2_RAW_FQ)."
        echo "  Stelle sicher, dass die Dekomprimierung erfolgreich war oder die .fq-Dateien dort vorhanden sind."
        echo "  Überspringe dieses Sample."
        continue # Springt zum nächsten Sample in der Schleife.
    fi

    # --- KORRIGIERTER CUTADAPT-AUFRUF ---
    # Hier werden Adapter entfernt und schlechte Qualität getrimmt.
    echo "  [2/5] Adapter und Quality Trimming mit cutadapt..."
    cutadapt \
        -a "$ADAPTER_SEQ" \
        -A "$ADAPTER_SEQ" \
        -q "$QUALITY_CUTOFF" \
        -m "$MIN_LENGTH_AFTER_TRIM" \
        --paired-output "$R2_TRIMMED_FQ" \
        -o "$R1_TRIMMED_FQ" \
        "$R1_RAW_FQ" "$R2_RAW_FQ" > "$CUTADAPT_LOG" 2>&1 # Die Umleitung > schreibt den Report in die Logdatei, 2>&1 auch Fehler.
    # **Bedeutung der Optionen (siehe PDF 5.1):**
    # -a: Adapter für Read1 (3'-Ende)
    # -A: Adapter für Read2 (3'-Ende)
    # -q: Qualitäts-Cutoff (hier 20)
    # -m: Minimale Länge nach Trimmen (hier 80)
    # --paired-output: Ausgabedatei für Read2 (wichtig für Paired-End)
    # -o: Ausgabedatei für Read1
    # "$R1_RAW_FQ" "$R2_RAW_FQ": Die beiden Eingabedateien.
    # **Deine Ausgabe zeigt:** Dieser Schritt hat jetzt funktioniert! Keine "command not found" Fehler mehr.

    # Sicherheitscheck: Wurden getrimmte Dateien erstellt und sind sie nicht leer?
    if [ ! -s "$R1_TRIMMED_FQ" ] || [ ! -s "$R2_TRIMMED_FQ" ]; then
        echo "  FEHLER: Cutadapt hat keine (oder leere) getrimmten Dateien für $SAMPLE_PREFIX erstellt."
        # ... (Fehlermeldung)
        continue
    fi
    # Deine Ausgabe zeigt diese Erfolgsmeldungen:
    echo "  Getrimmte Reads gespeichert in $R1_TRIMMED_FQ und $R2_TRIMMED_FQ"
    echo "  Cutadapt Log gespeichert in $CUTADAPT_LOG"

    # --- Mapping mit Hisat2 ---
    # Die getrimmten Reads werden jetzt auf das Referenzgenom gemappt.
    echo "  [3/5] Mapping mit Hisat2..."
    hisat2 \
        -p 4 \
        --dta \
        -x "$HISAT2_INDEX_PREFIX" \
        -1 "$R1_TRIMMED_FQ" \
        -2 "$R2_TRIMMED_FQ" \
        -S "$SAM_OUT" \
        --summary-file "$HISAT2_SUMMARY"
    # **Bedeutung der Optionen (siehe PDF 5.2):**
    # -p 4: Nutze 4 Prozessorkerne (Threads) – macht es schneller.
    # --dta: Option für spätere Transkriptom-Assemblierung (gute Praxis).
    # -x: Der Pfad zum Genom-Index (genome/Ahalleri).
    # -1, -2: Die getrimmten Input-Reads (R1 und R2).
    # -S: Die Ausgabe-SAM-Datei.
    # --summary-file: Eine Datei mit Mapping-Statistiken.
    # **Deine Ausgabe zeigt:** Für Noss.R1 z.B. "977311 reads; ... 74.22% overall alignment rate".
    # Das bedeutet, von den fast 1 Million Reads konnten ca. 74% auf das Genom gemappt werden. Das ist ein typischer Wert.

    # Sicherheitscheck SAM-Datei
    if [ ! -s "$SAM_OUT" ]; then
        # ... (Fehlermeldung)
        continue
    fi
    echo "  SAM-Datei gespeichert in $SAM_OUT"
    echo "  Hisat2 Summary gespeichert in $HISAT2_SUMMARY"

    # --- SAM zu BAM, Sortieren, Indexieren ---
    # SAM-Dateien sind groß. BAM ist komprimiert. Sortiert und indexiert wird für viele Tools gebraucht.
    echo "  [4/5] Konvertiere SAM zu BAM, sortiere und indexiere BAM..."
    samtools view -@ 4 -bS "$SAM_OUT" | samtools sort -@ 4 - -o "$SORTED_BAM_OUT"
    # samtools view -bS: Konvertiert SAM (-) zu BAM (-b). -@ 4 nutzt 4 Threads.
    # | (Pipe): Leitet die BAM-Ausgabe direkt an samtools sort.
    # samtools sort -o: Sortiert die BAM-Datei und speichert sie.
    samtools index "$SORTED_BAM_OUT" # Erstellt einen Index (.bai) für die BAM-Datei.
    # **Deine Ausgabe zeigt:** "[bam_sort_core] merging..." und Erfolgsmeldungen. Das ist gut.

    # SAM-Datei löschen, um Platz zu sparen.
    if [ -s "$SORTED_BAM_OUT" ]; then
        rm "$SAM_OUT"
        echo "  Sortierte BAM-Datei gespeichert in $SORTED_BAM_OUT und indexiert."
        echo "  Original SAM-Datei $SAM_OUT gelöscht."
    else
        # ... (Fehlermeldung)
        continue
    fi

    # --- Read Counting mit featureCounts ---
    # Zählt, wie viele Reads auf jedes Gen in der Annotationsdatei fallen.
    echo "  [5/5] Read Counting mit featureCounts..."
    featureCounts \
        -p \
        -t exon \
        -g gene_id \
        --primary \
        -a "$GENOME_ANNOTATION" \
        -o "${COUNT_FILE_PREFIX}.counts.tsv" \
        -T 4 \
        "$SORTED_BAM_OUT"
    # **Bedeutung der Optionen (siehe PDF 5.3):**
    # -p: Paired-End Daten.
    # -t exon: Zähle Reads, die auf "exon"-Features in der GFF-Datei mappen.
    # -g gene_id: Gruppiere die Exon-Counts nach dem Attribut "gene_id", um Gen-Counts zu erhalten.
    # --primary: Zähle nur primäre Alignments (die besten).
    # -a: Die Annotationsdatei (Ahalleri.gff).
    # -o: Die Ausgabedatei für die Zähltabelle.
    # -T 4: Nutze 4 Threads.
    # "$SORTED_BAM_OUT": Die sortierte BAM-Datei als Input.
    # **Deine Ausgabe zeigt:** Das featureCounts-Logo und Statistiken wie "Successfully assigned alignments : 736044 (34.7%)" für Noss.R1.
    # Das bedeutet, ca. 35% der gemappten Reads konnten eindeutig einem Gen zugeordnet werden. Das ist auch ein typischer Wert.

    # Summary-Datei umbenennen.
    if [ -f "${COUNT_FILE_PREFIX}.counts.tsv.summary" ]; then
        mv "${COUNT_FILE_PREFIX}.counts.tsv.summary" "${COUNT_FILE_PREFIX}.counts_summary.txt"
        echo "  Read Counts gespeichert in ${COUNT_FILE_PREFIX}.counts.tsv"
        echo "  FeatureCounts Summary gespeichert in ${COUNT_FILE_PREFIX}.counts_summary.txt"
    else
        echo "  WARNUNG: FeatureCounts Summary-Datei nicht gefunden..."
    fi

    echo "--- Sample $SAMPLE_PREFIX erfolgreich abgeschlossen ---"
done # Ende der Sample-Schleife

DURATION=$SECONDS
echo ""
echo "======================================================"
echo "=== RNA-SEQ PIPELINE ABGESCHLOSSEN ==="
echo "Gesamtdauer: $(($DURATION / 3600)) Stunden, $((($DURATION / 60) % 60)) Minuten und $(($DURATION % 60)) Sekunden"
echo "======================================================"
# **Deine Ausgabe zeigt:** Das Skript ist für alle 6 Samples durchgelaufen und hat ca. 3.5 Minuten gebraucht.
# **Hausaufgabe a) ist hiermit vollständig und erfolgreich erledigt!**
```

**Skript 2: `run_seqstats_analysis.sh` – Qualitätskontrolle**

**Was dieses Skript tut (Hausaufgabe b):**
Es verwendet das Werkzeug `seqstats`, um einfache Statistiken über deine Read-Dateien zu erstellen – einmal für alle getrimmten Dateien und als Bonus für eine rohe Datei zum Vergleich.

```bash
#!/bin/bash
# Shebang

echo "=== STARTING SEQSTATS ANALYSIS (Homework Version) ==="
SECONDS=0

# --- 0. Verzeichnisse und Pfade definieren ---
BASE_DIR=$(pwd)
TRIM_DIR="${BASE_DIR}/trimmed_reads"         # Wo die getrimmten Reads von Skript 1 liegen.
RAW_DIR_PROCESSED="${BASE_DIR}/raw_reads_all" # Wo die (dekomprimierten) Rohdaten liegen.
SEQSTATS_REPORTS_DIR="${BASE_DIR}/seqstats_reports" # Neuer Ordner für die Berichte.
mkdir -p "$SEQSTATS_REPORTS_DIR" # Erstellt den Ordner, falls nicht vorhanden.

# SAMPLES-Array, wie in Skript 1.
SAMPLES=("Noss.R1" "Noss.R2" "Noss.R3" "Pais.R1" "Pais.R2" "Pais.R3")

# --- 1. Erstelle Seqstats-Berichte für getrimmte Reads ---
echo "--- [1/2] Erstelle Seqstats-Berichte für getrimmte Reads ---"
# Schleife über alle Samples und dann über Read 1 und Read 2 jedes Samples.
for SAMPLE_PREFIX in "${SAMPLES[@]}"; do
    for PAIR_NUM in 1 2; do
        TRIMMED_FQ_FILE="$TRIM_DIR/${SAMPLE_PREFIX}.${PAIR_NUM}.trimmed.fq" # z.B. trimmed_reads/Noss.R1.1.trimmed.fq
        REPORT_OUT_TRIMMED="$SEQSTATS_REPORTS_DIR/${SAMPLE_PREFIX}.${PAIR_NUM}.trimmed.seqstats.txt" # z.B. seqstats_reports/Noss.R1.1.trimmed.seqstats.txt

        if [ -s "$TRIMMED_FQ_FILE" ]; then # Prüft, ob die getrimmte Datei existiert und nicht leer ist.
            echo "  Generiere Bericht für $TRIMMED_FQ_FILE -> $REPORT_OUT_TRIMMED"
            seqstats "$TRIMMED_FQ_FILE" > "$REPORT_OUT_TRIMMED" # Führt seqstats aus und leitet die Ausgabe in die Datei.
        else
            echo "  WARNUNG: Getrimmte Datei $TRIMMED_FQ_FILE nicht gefunden oder leer."
            # ...
        fi
    done
done
echo "--- Seqstats-Berichte für getrimmte Reads erstellt. ---"
# **Deine Ausgabe zeigt:** Für alle 12 getrimmten Dateien wurde ein Bericht generiert. Das ist korrekt.

# --- 2. Bonus: Seqstats für rohe Reads und Vergleich ---
echo ""
echo "--- [2/2] Bonus: Seqstats für rohe Reads (Beispiel Noss.R1.1) und Vergleich ---"
RAW_FQ_EXAMPLE_R1="$RAW_DIR_PROCESSED/Noss.R1.1.sub.fq" # Die rohe Beispieldatei
REPORT_RAW_OUT_EXAMPLE="$SEQSTATS_REPORTS_DIR/Noss.R1.1.raw.seqstats.txt" # Ihr Bericht
TRIMMED_FQ_EXAMPLE_R1="$TRIM_DIR/Noss.R1.1.trimmed.fq" # Die entsprechende getrimmte Datei
REPORT_TRIMMED_OUT_EXAMPLE="$SEQSTATS_REPORTS_DIR/Noss.R1.1.trimmed.seqstats.txt" # Ihr Bericht (schon erstellt)
MIN_LENGTH_CUTADAPT_PARAM=80 # Wichtig für die Diskussion!

if [ -s "$RAW_FQ_EXAMPLE_R1" ]; then
    echo "  Generiere Bericht für rohe Datei $RAW_FQ_EXAMPLE_R1 -> $REPORT_RAW_OUT_EXAMPLE"
    seqstats "$RAW_FQ_EXAMPLE_R1" > "$REPORT_RAW_OUT_EXAMPLE"

    echo ""
    echo "  Vergleich für Noss.R1.1 (Roh vs. Getrimmt):"
    echo "  ----------------------------------------------------"
    echo "  Rohe Daten ($RAW_FQ_EXAMPLE_R1):"
    if [ -s "$REPORT_RAW_OUT_EXAMPLE" ]; then
        cat "$REPORT_RAW_OUT_EXAMPLE" # cat gibt den Inhalt der Datei aus.
    else
        echo "    Bericht für rohe Daten ($REPORT_RAW_OUT_EXAMPLE) nicht gefunden."
    fi
    echo "  ----------------------------------------------------"
    echo "  Getrimmte Daten ($TRIMMED_FQ_EXAMPLE_R1):"
    if [ -s "$REPORT_TRIMMED_OUT_EXAMPLE" ]; then
        cat "$REPORT_TRIMMED_OUT_EXAMPLE"
    else
        echo "    Bericht für getrimmte Daten ($REPORT_TRIMMED_OUT_EXAMPLE) nicht gefunden."
    fi
    echo "  ----------------------------------------------------"
    echo ""
    # Der Diskussionsteil, der dir hilft, die Ergebnisse zu interpretieren:
    echo "  Diskussion (manuell ausfüllen basierend auf den obigen Werten):"
    # ... (Fragen und Erwartungen)
else
    echo "  WARNUNG: Rohe Beispieldatei $RAW_FQ_EXAMPLE_R1 nicht gefunden oder leer..."
fi

DURATION=$SECONDS
echo ""
echo "======================================================"
echo "=== SEQSTATS ANALYSIS ABGESCHLOSSEN ==="
echo "Gesamtdauer: $(($DURATION / 60)) Minuten und $(($DURATION % 60)) Sekunden"
echo "======================================================"
# **Deine Ausgabe zeigt:**
# Rohe Daten Noss.R1.1: 1 Mio Reads, alle 150bp lang.
# Getrimmte Daten Noss.R1.1: 977.311 Reads (weniger!), durchschnittlich 138.32bp, min. Länge 80bp.
# **Das ist genau das, was wir erwarten!** Das Trimmen hat Adapter/schlechte Qualität entfernt und Reads,
# die kürzer als 80bp wurden, verworfen. Die minimale Länge von 80bp bei den getrimmten Reads
# bestätigt, dass der `-m 80` Parameter von cutadapt korrekt angewendet wurde.
# **Hausaufgabe b) ist hiermit vollständig und erfolgreich erledigt!**
```

**Skript 3: `create_submission_zip.sh` – Ergebnisse verpacken**

**Was dieses Skript tut (Hausaufgabe c):**
Es sammelt die 12 (getrimmten) `seqstats`-Berichte und die 6 `featureCounts`-Zähltabellen und packt sie in eine ZIP-Datei.

```bash
#!/bin/bash
# Shebang

echo "=== CREATING SUBMISSION ZIP ARCHIVE ==="
# ... (Pfade werden definiert) ...
SUBMISSION_ZIP_FILE="${BASE_DIR}/rna_seq_results_$(whoami)_$(date +%Y%m%d).zip"
# Der Dateiname des ZIPs enthält deinen Benutzernamen und das Datum.

# ... (Prüfungen, ob Ordner existieren) ...

echo "Folgende Dateien werden zum ZIP-Archiv hinzugefügt:"
FILES_TO_ZIP=() # Eine leere Liste für die Dateinamen.

echo "Suche Seqstats-Berichte (*.seqstats.txt) in $SEQSTATS_REPORTS_DIR..."
report_files=("$SEQSTATS_REPORTS_DIR"/*.seqstats.txt) # Findet alle .seqstats.txt Dateien.
if [ ${#report_files[@]} -gt 0 ] && [ -e "${report_files[0]}" ]; then
    for file in "${report_files[@]}"; do
        echo "  - $(basename "$SEQSTATS_REPORTS_DIR")/$(basename "$file")" # Zeigt den relativen Pfad an.
        FILES_TO_ZIP+=("$(basename "$SEQSTATS_REPORTS_DIR")/$(basename "$file")") # Fügt zur Liste hinzu.
    done
else
    echo "WARNUNG: Keine Seqstats-Berichte (*.seqstats.txt) in $SEQSTATS_REPORTS_DIR gefunden!"
fi

echo "Suche Count-Dateien (*.counts.tsv) in $COUNT_DIR..."
count_files=("$COUNT_DIR"/*.counts.tsv) # Findet alle .counts.tsv Dateien.
if [ ${#count_files[@]} -gt 0 ] && [ -e "${count_files[0]}" ]; then
    for file in "${count_files[@]}"; do
        echo "  - $(basename "$COUNT_DIR")/$(basename "$file")"
        FILES_TO_ZIP+=("$(basename "$COUNT_DIR")/$(basename "$file")")
    done
else
    echo "WARNUNG: Keine Count-Dateien (*.counts.tsv) in $COUNT_DIR gefunden!"
fi

if [ ${#FILES_TO_ZIP[@]} -eq 0 ]; then
    echo "FEHLER: Keine Dateien zum Zippen gefunden. Breche ab."
    exit 1
fi

echo "Erstelle ZIP-Archiv: $SUBMISSION_ZIP_FILE"
( # Start einer Subshell, damit cd das äußere Skript nicht beeinflusst
  cd "$BASE_DIR" || exit 1 # Wechselt ins Basisverzeichnis, damit Pfade im ZIP relativ sind.
  zip "$SUBMISSION_ZIP_FILE" "${FILES_TO_ZIP[@]}" # Der eigentliche zip-Befehl.
)

if [ -f "$SUBMISSION_ZIP_FILE" ]; then
    echo "ZIP-Archiv erfolgreich erstellt: $SUBMISSION_ZIP_FILE"
    echo "Inhalt des ZIP-Archivs:"
    unzip -l "$SUBMISSION_ZIP_FILE" # Zeigt den Inhalt des ZIPs.
else
    echo "FEHLER: ZIP-Archiv konnte nicht erstellt werden."
    exit 1
fi

echo "=== ZIPPING ABGESCHLOSSEN ==="
# **Deine Ausgabe zeigt:**
# Es wurden 13 seqstats-Berichte gefunden (inklusive des Noss.R1.1.raw.seqstats.txt) und 6 count-Dateien.
# Das ZIP-Archiv wurde erfolgreich erstellt und enthält diese 19 Dateien.
# Die Pfade im ZIP sind relativ (z.B. seqstats_reports/Noss.R1.1.raw.seqstats.txt), was gut ist.
# **Hausaufgabe c) ist hiermit vollständig und erfolgreich erledigt!**
```

---

**Teil 3: Hausaufgabe d) – Der Text**

Du wurdest gebeten, die wichtigsten Schritte der RNA-Seq-Analyse in eigenen Worten zu erklären. Hier ist eine Erinnerung an die Kernschritte, die du in deinem Text erwähnen solltest, basierend auf dem, was die Skripte tun und was in den PDFs/Folien steht:

1.  **(Optional, aber gute Praxis) Qualitätskontrolle der Rohdaten:** Bevor man anfängt, schaut man sich die Qualität der ursprünglichen Sequenzierdaten an. (Werkzeuge: FastQC, MultiQC, oder wie hier `seqstats` für eine einfache Übersicht).
2.  **Adapter-Trimming und Qualitätsfilterung:** Entfernen von technischen Adaptersequenzen und Basen niedriger Qualität von den Reads. Zu kurze Reads werden verworfen. (Werkzeug: `cutadapt`).
3.  **Mapping/Alignment:** Die sauberen Reads werden an ein Referenzgenom (oder Transkriptom) ausgerichtet, um ihre Herkunft zu bestimmen. (Werkzeug: `Hisat2`).
4.  **Quantifizierung (Read Counting):** Zählen, wie viele Reads jedem Gen zugeordnet wurden. Dies ist ein Maß für die Genexpression. (Werkzeug: `featureCounts`).
5.  **(Nächster Schritt, nicht in diesen Skripten) Differentielle Genexpressionsanalyse:** Vergleichen der Genexpressionslevel zwischen verschiedenen Bedingungen (z.B. Noss vs. Pais), um signifikante Unterschiede zu finden.

---

**Zusammenfassung und Fazit:**

Du hast erfolgreich eine komplette, wenn auch grundlegende, RNA-Seq-Analyse-Pipeline für sechs Probenpaare automatisiert!

*   **Datenvorbereitung:** Deine Rohdaten (ob komprimiert oder schon dekomprimiert in `raw_reads_all`) wurden korrekt vom ersten Skript erkannt und für die Analyse bereitgestellt.
*   **Trimming:** `cutadapt` hat mit den wichtigen Parametern aus der PDF (Adapter, Qualitätsschnitt 20, **Minimallänge 80 Basen**) funktioniert. Der `seqstats`-Vergleich hat dies bestätigt.
*   **Mapping:** `Hisat2` hat einen guten Anteil deiner Reads auf das `Ahalleri`-Genom gemappt.
*   **Counting:** `featureCounts` hat die Reads den Genen zugeordnet und Zähltabellen erstellt.
*   **Qualitätskontrolle:** `seqstats` hat dir geholfen, die Auswirkungen des Trimmings zu sehen.
*   **Organisation:** Alle Ausgabedateien wurden in den richtigen Ordnern gespeichert, und die für die Abgabe relevanten Dateien wurden korrekt gezippt.

**Deine Hausaufgaben sind damit im Wesentlichen erledigt!** jetzt nur noch:
1.  Deine Beobachtungen zum `seqstats`-Vergleich (Hausaufgabe b, Diskussionsteil) notieren.
2.  Den Text für Hausaufgabe d) verfassen.
3.  Das ZIP-Archiv und den Text per E-Mail einreichen.


**Hausaufgabe b) – Deine Beobachtungen zum `seqstats`-Vergleich (Noss.R1.1)**

Hier ist ein Entwurf für deine Beobachtungen, der direkt auf den Zahlen aus deiner Terminalausgabe basiert:

```
  Deine Beobachtungen:
  - **Anzahl Sequenzen (Total n):** Die Anzahl der Sequenzen ist von 1.000.000 in den rohen Daten auf 977.311 in den getrimmten Daten gesunken. Dies entspricht den Erwartungen, da cutadapt Reads (und ihre Paare) verwirft, die nach dem Trimmen kürzer als die vorgegebene minimale Länge von 80 Basen sind. Etwa 2,27% der Reads wurden also verworfen.

  - **Gesamtzahl Basen (Total seq):** Die Gesamtzahl der Basen hat sich von 150.000.000 bp (Rohdaten) auf 135.182.567 bp (getrimmte Daten) reduziert. Dies ist ebenfalls erwartet, da sowohl Adaptersequenzen als auch Basen niedriger Qualität von den Enden der Reads entfernt wurden, was die Gesamtmenge an Sequenzinformation verringert.

  - **Durchschnittliche Sequenzlänge (Avg. seq):** Die durchschnittliche Sequenzlänge ist von 150.00 bp auf 138.32 bp gesunken. Dies ist eine direkte Folge des Trimmings von Adaptern und Low-Quality-Basen.

  - **Minimale / Maximale Sequenzlänge (Min seq / Max seq):**
    - In den Rohdaten war die minimale und maximale Länge einheitlich 150 bp.
    - In den getrimmten Daten ist die **minimale Länge jetzt 80 bp**. Dies bestätigt exakt die Anwendung des Parameters `-m 80` in cutadapt. Reads, die kürzer geworden wären, wurden verworfen. Die maximale Länge ist bei 150 bp geblieben, was bedeutet, dass einige Reads entweder gar nicht oder nur minimal getrimmt werden mussten und ihre ursprüngliche Länge (oder eine Länge nahe daran) beibehalten haben.

  - **Allgemeiner Effekt:** Das Trimmen und Clippen hatte einen deutlichen und erwarteten Effekt. Es hat die Datenmenge reduziert, potenziell problematische Sequenzabschnitte (Adapter, niedrige Qualität) entfernt und sichergestellt, dass alle verbleibenden Reads eine Mindestlänge für die weitere Analyse aufweisen. Die Ergebnisse entsprechen voll und ganz den Erwartungen basierend auf den verwendeten cutadapt-Parametern.
```

---

**Hausaufgabe d) – Text zur RNA-Seq Analyse**

Hier ist ein Entwurf für deinen Text. Er basiert auf den Schritten, die deine Skripte durchgeführt haben, und den Werkzeugen, die du verwendet hast. Formuliere ihn gerne noch etwas um, damit er wirklich "deine eigenen Worte" widerspiegelt, aber die Kerninformationen sind enthalten.

```

**Zusammenfassung der wichtigsten Schritte der RNA-Seq Analyse**

Die Analyse von RNA-Sequenzierungsdaten (RNA-Seq) ist ein mehrstufiger Prozess, der darauf abzielt, die Genexpression in biologischen Proben zu quantifizieren und zu vergleichen. Die von mir durchgeführte grundlegende Pipeline umfasste die folgenden Hauptschritte:

1.  **Datenvorbereitung (Dekompression):**
    Die Rohdaten vom Sequenzierer liegen oft in komprimierter Form vor (z.B. als `.fq.gz` oder `.fq.zip` Dateien). Der erste Schritt bestand darin, diese Dateien zu dekomprimieren, um unkomprimierte FASTQ-Dateien (`.fq`) für die weitere Verarbeitung zu erhalten. Dies wurde im Skript durch `gunzip` und `unzip` realisiert.

2.  **Adapter-Trimming und Qualitätsfilterung:**
    Sequenzier-Reads können technische Adaptersequenzen enthalten, die während der Bibliothekspräparation hinzugefügt wurden, sowie Basen von niedriger Qualität, oft an den 3'-Enden der Reads. Diese müssen entfernt werden, um die Genauigkeit der nachfolgenden Analyseschritte zu verbessern.
    *   **Werkzeug:** Für diese Aufgabe wurde `cutadapt` verwendet.
    *   **Parameter:** Wichtige Parameter waren die spezifische Adaptersequenz (`AGATCGGAAGAG`), ein Qualitäts-Cutoff von `20` zum Entfernen von Basen niedriger Qualität vom 3'-Ende und eine minimale Leselänge von `80` Basen nach dem Trimmen, um sicherzustellen, dass nur ausreichend lange Reads für das Mapping verwendet werden. Die Verarbeitung erfolgte im Paired-End-Modus.

3.  **Mapping (Alignment) der Reads auf ein Referenzgenom:**
    Die bereinigten und getrimmten Reads wurden dann gegen ein Referenzgenom (in diesem Fall *Arabidopsis halleri*) alignt. Ziel dieses Schrittes ist es, die ursprüngliche genomische Position jedes Reads zu bestimmen.
    *   **Werkzeug:** Für das Mapping wurde `Hisat2` eingesetzt.
    *   **Voraussetzung:** `Hisat2` benötigt einen vorab erstellten Index des Referenzgenoms, um das Alignment effizient durchführen zu können.
    *   **Ausgabe:** Das Ergebnis des Mappings wurde zunächst als SAM-Datei (Sequence Alignment/Map) gespeichert.

4.  **Konvertierung und Aufbereitung der Mapping-Dateien:**
    SAM-Dateien sind textbasiert und können sehr groß sein. Daher wurden sie in das komprimierte Binärformat BAM (Binary Alignment/Map) konvertiert. Anschließend wurden die BAM-Dateien nach genomischen Koordinaten sortiert und indexiert. Diese Schritte sind oft notwendig für nachfolgende Analysen.
    *   **Werkzeug:** Für diese Aufgaben wurde `samtools` (mit den Befehlen `view`, `sort` und `index`) verwendet.

5.  **Quantifizierung der Genexpression (Read Counting):**
    Im letzten Schritt der Pipeline wurde gezählt, wie viele der gemappten Reads jedem annotierten Gen im Referenzgenom zugeordnet werden können. Die Anzahl der Reads, die einem Gen zugeordnet werden, dient als Maß für dessen Expressionslevel in der jeweiligen Probe.
    *   **Werkzeug:** Für das Read Counting wurde `featureCounts` (aus dem Subread-Paket) verwendet.
    *   **Voraussetzung:** `featureCounts` benötigt eine Genom-Annotationsdatei (hier im GFF-Format: `Ahalleri.gff`), die die Positionen und Strukturen der Gene definiert.
    *   **Ausgabe:** Das Ergebnis ist eine Zähltabelle (Counts pro Gen pro Sample), die die Grundlage für weiterführende Analysen wie die differentielle Genexpression bildet.

Zusätzlich wurde mit dem Werkzeug `seqstats` eine Qualitätskontrolle der Roh- und getrimmten Daten durchgeführt, um die Effektivität des Trimming-Schritts zu bewerten. Alle Schritte wurden für sechs Paired-End-Probenpaare automatisiert in einer Schleife durchgeführt.

