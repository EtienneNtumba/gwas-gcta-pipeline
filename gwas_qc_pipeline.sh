#!/bin/bash

################################################################################
# Pipeline QC GWAS Robuste et Efficace
# Pour: 3210 échantillons
# Objectif: GWAS avec EMMAX + Calcul PRS de qualité
# Auteur: Pipeline optimisé pour données génomiques
################################################################################

# CONFIGURATION
THREADS=8
INPUT_PREFIX="/mnt/lustre/groups/CBBI0818/TZ/data/TZ.HbF.phased.3"  # Modifier selon vos fichiers PLINK (.bed/.bim/.fam)
OUTPUT_DIR="/mnt/lustre/groups/CBBI0818/TZ/data/qc_results"
FINAL_PREFIX="qc_final"

# Seuils de QC (ajustables selon votre étude)
MAF_THRESHOLD=0.01          # MAF minimum (1% pour GWAS)
GENO_SNP=0.02              # Taux de missing maximum par SNP (2%)
GENO_INDIV=0.02            # Taux de missing maximum par individu (2%)
HWE_PVALUE=1e-6            # p-value HWE (stringent)
HET_SD=3                   # Écart-type pour hétérozygotie (±3 SD)
IBD_THRESHOLD=0.1875       # Seuil pour identifier relatives (> 2nd degré)
LD_R2=0.2                  # R² pour pruning LD
LD_WINDOW=50               # Fenêtre pour LD (kb)
LD_STEP=5                  # Step pour LD

# Créer les dossiers
mkdir -p ${OUTPUT_DIR}/{step1_initial,step2_snp_qc,step3_sample_qc,step4_population,step5_final,logs,plots}

echo "=================================="
echo "GWAS QC Pipeline - Démarrage"
echo "Échantillons: 3210"
echo "Date: $(date)"
echo "=================================="

################################################################################
# ÉTAPE 0: VÉRIFICATION INITIALE DES DONNÉES
################################################################################
echo "[ÉTAPE 0] Vérification initiale des données..."

plink --bfile ${INPUT_PREFIX} \
      --missing \
      --freq \
      --hardy \
      --het \
      --out ${OUTPUT_DIR}/step1_initial/initial_stats \
      --threads ${THREADS}

# Statistiques de base
plink --bfile ${INPUT_PREFIX} \
      --check-sex \
      --out ${OUTPUT_DIR}/step1_initial/sex_check \
      --threads ${THREADS}

echo "Statistiques initiales générées."

################################################################################
# ÉTAPE 1: QC DES VARIANTS (SNPs)
################################################################################
echo ""
echo "[ÉTAPE 1] Contrôle qualité des SNPs..."

# 1.1 - Filtrer sur le taux de génotypage par SNP
echo "  1.1 - Filtrage call rate SNPs (>${GENO_SNP})..."
plink --bfile ${INPUT_PREFIX} \
      --geno ${GENO_SNP} \
      --make-bed \
      --out ${OUTPUT_DIR}/step2_snp_qc/qc1_geno \
      --threads ${THREADS}

# 1.2 - Filtrer sur MAF
echo "  1.2 - Filtrage MAF (>${MAF_THRESHOLD})..."
plink --bfile ${OUTPUT_DIR}/step2_snp_qc/qc1_geno \
      --maf ${MAF_THRESHOLD} \
      --make-bed \
      --out ${OUTPUT_DIR}/step2_snp_qc/qc2_maf \
      --threads ${THREADS}

# 1.3 - Filtrer sur équilibre Hardy-Weinberg (HWE)
echo "  1.3 - Filtrage HWE (p>${HWE_PVALUE})..."
plink --bfile ${OUTPUT_DIR}/step2_snp_qc/qc2_maf \
      --hwe ${HWE_PVALUE} \
      --make-bed \
      --out ${OUTPUT_DIR}/step2_snp_qc/qc3_hwe \
      --threads ${THREADS}

echo "QC SNPs terminé."

################################################################################
# ÉTAPE 2: QC DES ÉCHANTILLONS (SAMPLES)
################################################################################
echo ""
echo "[ÉTAPE 2] Contrôle qualité des échantillons..."

# 2.1 - Filtrer sur le taux de génotypage par individu
echo "  2.1 - Filtrage call rate individus (>${GENO_INDIV})..."
plink --bfile ${OUTPUT_DIR}/step2_snp_qc/qc3_hwe \
      --mind ${GENO_INDIV} \
      --make-bed \
      --out ${OUTPUT_DIR}/step3_sample_qc/qc1_mind \
      --threads ${THREADS}

# 2.2 - Vérification du sexe (identifier discordances)
echo "  2.2 - Vérification du sexe..."
plink --bfile ${OUTPUT_DIR}/step3_sample_qc/qc1_mind \
      --check-sex \
      --out ${OUTPUT_DIR}/step3_sample_qc/sex_check \
      --threads ${THREADS}

# Créer liste des individus avec problèmes de sexe (PROBLEM dans le fichier .sexcheck)
awk '$5=="PROBLEM" {print $1, $2}' ${OUTPUT_DIR}/step3_sample_qc/sex_check.sexcheck > \
    ${OUTPUT_DIR}/step3_sample_qc/sex_discordance.txt

# Retirer les individus avec discordance de sexe
if [ -s ${OUTPUT_DIR}/step3_sample_qc/sex_discordance.txt ]; then
    echo "  Retrait de $(wc -l < ${OUTPUT_DIR}/step3_sample_qc/sex_discordance.txt) individus avec discordance de sexe"
    plink --bfile ${OUTPUT_DIR}/step3_sample_qc/qc1_mind \
          --remove ${OUTPUT_DIR}/step3_sample_qc/sex_discordance.txt \
          --make-bed \
          --out ${OUTPUT_DIR}/step3_sample_qc/qc2_sex \
          --threads ${THREADS}
else
    echo "  Aucune discordance de sexe détectée"
    cp ${OUTPUT_DIR}/step3_sample_qc/qc1_mind.bed ${OUTPUT_DIR}/step3_sample_qc/qc2_sex.bed
    cp ${OUTPUT_DIR}/step3_sample_qc/qc1_mind.bim ${OUTPUT_DIR}/step3_sample_qc/qc2_sex.bim
    cp ${OUTPUT_DIR}/step3_sample_qc/qc1_mind.fam ${OUTPUT_DIR}/step3_sample_qc/qc2_sex.fam
fi

# 2.3 - Hétérozygotie (identifier outliers)
echo "  2.3 - Calcul de l'hétérozygotie..."

# D'abord, faire un pruning LD pour hétérozygotie
plink --bfile ${OUTPUT_DIR}/step3_sample_qc/qc2_sex \
      --indep-pairwise ${LD_WINDOW} ${LD_STEP} ${LD_R2} \
      --out ${OUTPUT_DIR}/step3_sample_qc/ld_pruned \
      --threads ${THREADS}

# Calculer hétérozygotie sur SNPs indépendants
plink --bfile ${OUTPUT_DIR}/step3_sample_qc/qc2_sex \
      --extract ${OUTPUT_DIR}/step3_sample_qc/ld_pruned.prune.in \
      --het \
      --out ${OUTPUT_DIR}/step3_sample_qc/het \
      --threads ${THREADS}

# Créer script R pour identifier outliers d'hétérozygotie
cat > ${OUTPUT_DIR}/step3_sample_qc/het_outliers.R << 'RSCRIPT'
#!/usr/bin/env Rscript
het <- read.table("/mnt/lustre/groups/CBBI0818/TZ/data/qc_results/step3_sample_qc/het.het", header=TRUE)
het$HET_RATE <- (het$N.NM. - het$O.HOM.) / het$N.NM.
mean_het <- mean(het$HET_RATE)
sd_het <- sd(het$HET_RATE)

# Identifier outliers (±3 SD)
het$OUTLIER <- abs(het$HET_RATE - mean_het) > 3 * sd_het
outliers <- het[het$OUTLIER == TRUE, c("FID", "IID")]

write.table(outliers, "/mnt/lustre/groups/CBBI0818/TZ/data/qc_results/step3_sample_qc/het_outliers.txt", 
            quote=FALSE, row.names=FALSE, col.names=FALSE)

# Créer plot
png("/mnt/lustre/groups/CBBI0818/TZ/data/qc_results/plots/heterozygosity_distribution.png", width=800, height=600)
hist(het$HET_RATE, breaks=50, main="Distribution de l'hétérozygotie", 
     xlab="Taux d'hétérozygotie", col="lightblue")
abline(v=mean_het + 3*sd_het, col="red", lty=2, lwd=2)
abline(v=mean_het - 3*sd_het, col="red", lty=2, lwd=2)
abline(v=mean_het, col="blue", lwd=2)
dev.off()

cat("Moyenne hétérozygotie:", mean_het, "\n")
cat("SD:", sd_het, "\n")
cat("Outliers:", nrow(outliers), "\n")
RSCRIPT

chmod +x ${OUTPUT_DIR}/step3_sample_qc/het_outliers.R
Rscript ${OUTPUT_DIR}/step3_sample_qc/het_outliers.R

# Retirer outliers d'hétérozygotie
if [ -s ${OUTPUT_DIR}/step3_sample_qc/het_outliers.txt ]; then
    echo "  Retrait de $(wc -l < ${OUTPUT_DIR}/step3_sample_qc/het_outliers.txt) outliers d'hétérozygotie"
    plink --bfile ${OUTPUT_DIR}/step3_sample_qc/qc2_sex \
          --remove ${OUTPUT_DIR}/step3_sample_qc/het_outliers.txt \
          --make-bed \
          --out ${OUTPUT_DIR}/step3_sample_qc/qc3_het \
          --threads ${THREADS}
else
    echo "  Aucun outlier d'hétérozygotie"
    cp ${OUTPUT_DIR}/step3_sample_qc/qc2_sex.bed ${OUTPUT_DIR}/step3_sample_qc/qc3_het.bed
    cp ${OUTPUT_DIR}/step3_sample_qc/qc2_sex.bim ${OUTPUT_DIR}/step3_sample_qc/qc3_het.bim
    cp ${OUTPUT_DIR}/step3_sample_qc/qc2_sex.fam ${OUTPUT_DIR}/step3_sample_qc/qc3_het.fam
fi

echo "QC échantillons terminé."

################################################################################
# ÉTAPE 3: DÉTECTION DE PARENTÉ ET DUPLICATS (IBD/IBS)
################################################################################
echo ""
echo "[ÉTAPE 3] Détection de parenté et duplicats..."

# 3.1 - Pruning LD pour calcul IBD
echo "  3.1 - LD pruning pour IBD..."
plink --bfile ${OUTPUT_DIR}/step3_sample_qc/qc3_het \
      --indep-pairwise ${LD_WINDOW} ${LD_STEP} ${LD_R2} \
      --out ${OUTPUT_DIR}/step4_population/ld_pruned_ibd \
      --threads ${THREADS}

# 3.2 - Calcul IBD (Identity By Descent)
echo "  3.2 - Calcul IBD..."
plink --bfile ${OUTPUT_DIR}/step3_sample_qc/qc3_het \
      --extract ${OUTPUT_DIR}/step4_population/ld_pruned_ibd.prune.in \
      --genome \
      --min 0.1875 \
      --out ${OUTPUT_DIR}/step4_population/ibd \
      --threads ${THREADS}

# 3.3 - Identifier paires de relatives
# Créer script pour gérer les relatives
cat > ${OUTPUT_DIR}/step4_population/identify_relatives.R << 'RSCRIPT'
#!/usr/bin/env Rscript
if (file.exists("/mnt/lustre/groups/CBBI0818/TZ/data/qc_results/step4_population/ibd.genome")) {
  ibd <- read.table("/mnt/lustre/groups/CBBI0818/TZ/data/qc_results/step4_population/ibd.genome", header=TRUE)
  
  # Relatives: PI_HAT > 0.1875 (2nd degré ou plus proche)
  relatives <- ibd[ibd$PI_HAT > 0.1875, ]
  
  if (nrow(relatives) > 0) {
    # Pour chaque paire, garder l'individu avec le meilleur call rate
    # On va simplement prendre le premier de chaque paire (IID1)
    to_remove <- unique(relatives$IID1)
    to_remove_df <- data.frame(FID=relatives$FID1[match(to_remove, relatives$IID1)], 
                                IID=to_remove)
    
    write.table(to_remove_df, "/mnt/lustre/groups/CBBI0818/TZ/data/qc_results/step4_population/relatives_to_remove.txt",
                quote=FALSE, row.names=FALSE, col.names=FALSE)
    cat("Paires de relatives détectées:", nrow(relatives), "\n")
    cat("Individus à retirer:", nrow(to_remove_df), "\n")
  } else {
    cat("Aucune paire de relatives détectée\n")
    file.create("/mnt/lustre/groups/CBBI0818/TZ/data/qc_results/step4_population/relatives_to_remove.txt")
  }
} else {
  cat("Aucun fichier IBD (aucune paire de relatives)\n")
  file.create("/mnt/lustre/groups/CBBI0818/TZ/data/qc_results/step4_population/relatives_to_remove.txt")
}
RSCRIPT

chmod +x ${OUTPUT_DIR}/step4_population/identify_relatives.R
Rscript ${OUTPUT_DIR}/step4_population/identify_relatives.R

# Retirer relatives
if [ -s ${OUTPUT_DIR}/step4_population/relatives_to_remove.txt ]; then
    echo "  Retrait de $(wc -l < ${OUTPUT_DIR}/step4_population/relatives_to_remove.txt) individus apparentés"
    plink --bfile ${OUTPUT_DIR}/step3_sample_qc/qc3_het \
          --remove ${OUTPUT_DIR}/step4_population/relatives_to_remove.txt \
          --make-bed \
          --out ${OUTPUT_DIR}/step4_population/qc4_ibd \
          --threads ${THREADS}
else
    echo "  Aucun individu apparenté à retirer"
    cp ${OUTPUT_DIR}/step3_sample_qc/qc3_het.bed ${OUTPUT_DIR}/step4_population/qc4_ibd.bed
    cp ${OUTPUT_DIR}/step3_sample_qc/qc3_het.bim ${OUTPUT_DIR}/step4_population/qc4_ibd.bim
    cp ${OUTPUT_DIR}/step3_sample_qc/qc3_het.fam ${OUTPUT_DIR}/step4_population/qc4_ibd.fam
fi

################################################################################
# ÉTAPE 4: ANALYSE DE STRUCTURE DE POPULATION (PCA)
################################################################################
echo ""
echo "[ÉTAPE 4] Analyse de la structure de population (PCA)..."

# 4.1 - LD pruning pour PCA
echo "  4.1 - LD pruning pour PCA..."
plink --bfile ${OUTPUT_DIR}/step4_population/qc4_ibd \
      --indep-pairwise ${LD_WINDOW} ${LD_STEP} ${LD_R2} \
      --out ${OUTPUT_DIR}/step4_population/ld_pruned_pca \
      --threads ${THREADS}

# 4.2 - Calcul PCA
echo "  4.2 - Calcul des composantes principales..."
plink --bfile ${OUTPUT_DIR}/step4_population/qc4_ibd \
      --extract ${OUTPUT_DIR}/step4_population/ld_pruned_pca.prune.in \
      --pca 20 \
      --out ${OUTPUT_DIR}/step4_population/pca \
      --threads ${THREADS}

# 4.3 - Visualisation PCA
cat > ${OUTPUT_DIR}/step4_population/plot_pca.R << 'RSCRIPT'
#!/usr/bin/env Rscript
pca <- read.table("/mnt/lustre/groups/CBBI0818/TZ/data/qc_results/step4_population/pca.eigenvec", header=FALSE)
eigenval <- read.table("/mnt/lustre/groups/CBBI0818/TZ/data/qc_results/step4_population/pca.eigenval", header=FALSE)

# Calculer variance expliquée
var_explained <- eigenval$V1 / sum(eigenval$V1) * 100

# Plot PC1 vs PC2
png("/mnt/lustre/groups/CBBI0818/TZ/data/qc_results/plots/pca_pc1_pc2.png", width=1000, height=800)
plot(pca$V3, pca$V4, 
     xlab=paste0("PC1 (", round(var_explained[1], 2), "%)"),
     ylab=paste0("PC2 (", round(var_explained[2], 2), "%)"),
     main="Analyse en Composantes Principales (PC1 vs PC2)",
     pch=19, col=rgb(0,0,1,0.5))
dev.off()

# Plot PC1 vs PC3
png("/mnt/lustre/groups/CBBI0818/TZ/data/qc_results/plots/pca_pc1_pc3.png", width=1000, height=800)
plot(pca$V3, pca$V5, 
     xlab=paste0("PC1 (", round(var_explained[1], 2), "%)"),
     ylab=paste0("PC3 (", round(var_explained[3], 2), "%)"),
     main="Analyse en Composantes Principales (PC1 vs PC3)",
     pch=19, col=rgb(0,0,1,0.5))
dev.off()

# Scree plot
png("/mnt/lustre/groups/CBBI0818/TZ/data/qc_results/plots/pca_scree.png", width=1000, height=800)
barplot(var_explained[1:10], 
        names.arg=1:10,
        xlab="Composante Principale",
        ylab="Variance Expliquée (%)",
        main="Scree Plot - Variance Expliquée par PC",
        col="steelblue")
dev.off()

cat("\nVariance expliquée par les 10 premières PCs:\n")
print(data.frame(PC=1:10, Variance=round(var_explained[1:10], 2)))
RSCRIPT

chmod +x ${OUTPUT_DIR}/step4_population/plot_pca.R
Rscript ${OUTPUT_DIR}/step4_population/plot_pca.R

# 4.4 - Détection d'outliers de population (optionnel mais recommandé)
echo "  4.3 - Détection d'outliers de population..."
cat > ${OUTPUT_DIR}/step4_population/detect_pop_outliers.R << 'RSCRIPT'
#!/usr/bin/env Rscript
pca <- read.table("/mnt/lustre/groups/CBBI0818/TZ/data/qc_results/step4_population/pca.eigenvec", header=FALSE)

# Utiliser les 4 premières PCs pour détecter outliers
pc_data <- pca[, 3:6]

# Outliers = individus à plus de 6 SD sur n'importe quelle PC
outliers <- c()
for (i in 1:4) {
  mean_pc <- mean(pc_data[,i])
  sd_pc <- sd(pc_data[,i])
  outlier_idx <- which(abs(pc_data[,i] - mean_pc) > 6 * sd_pc)
  outliers <- c(outliers, outlier_idx)
}

outliers <- unique(outliers)

if (length(outliers) > 0) {
  outliers_df <- pca[outliers, 1:2]
  write.table(outliers_df, "/mnt/lustre/groups/CBBI0818/TZ/data/qc_results/step4_population/pop_outliers.txt",
              quote=FALSE, row.names=FALSE, col.names=FALSE)
  cat("Outliers de population détectés:", length(outliers), "\n")
} else {
  cat("Aucun outlier de population détecté\n")
  file.create("/mnt/lustre/groups/CBBI0818/TZ/data/qc_results/step4_population/pop_outliers.txt")
}
RSCRIPT

chmod +x ${OUTPUT_DIR}/step4_population/detect_pop_outliers.R
Rscript ${OUTPUT_DIR}/step4_population/detect_pop_outliers.R

# Retirer outliers de population si détectés
if [ -s ${OUTPUT_DIR}/step4_population/pop_outliers.txt ]; then
    echo "  Retrait de $(wc -l < ${OUTPUT_DIR}/step4_population/pop_outliers.txt) outliers de population"
    plink --bfile ${OUTPUT_DIR}/step4_population/qc4_ibd \
          --remove ${OUTPUT_DIR}/step4_population/pop_outliers.txt \
          --make-bed \
          --out ${OUTPUT_DIR}/step5_final/${FINAL_PREFIX} \
          --threads ${THREADS}
else
    echo "  Aucun outlier de population à retirer"
    cp ${OUTPUT_DIR}/step4_population/qc4_ibd.bed ${OUTPUT_DIR}/step5_final/${FINAL_PREFIX}.bed
    cp ${OUTPUT_DIR}/step4_population/qc4_ibd.bim ${OUTPUT_DIR}/step5_final/${FINAL_PREFIX}.bim
    cp ${OUTPUT_DIR}/step4_population/qc4_ibd.fam ${OUTPUT_DIR}/step5_final/${FINAL_PREFIX}.fam
fi

################################################################################
# ÉTAPE 5: DONNÉES FINALES ET STATISTIQUES
################################################################################
echo ""
echo "[ÉTAPE 5] Génération des données finales..."

# 5.1 - Statistiques finales
plink --bfile ${OUTPUT_DIR}/step5_final/${FINAL_PREFIX} \
      --missing \
      --freq \
      --hardy \
      --het \
      --out ${OUTPUT_DIR}/step5_final/final_stats \
      --threads ${THREADS}

# 5.2 - Créer fichier de covariables (sexe + 10 PCs)
echo "  5.1 - Création du fichier de covariables..."
cat > ${OUTPUT_DIR}/step5_final/create_covariates.R << 'RSCRIPT'
#!/usr/bin/env Rscript
# Charger PCA
pca <- read.table("/mnt/lustre/groups/CBBI0818/TZ/data/qc_results/step4_population/pca.eigenvec", header=FALSE)
colnames(pca) <- c("FID", "IID", paste0("PC", 1:20))

# Charger info de sexe
fam <- read.table("/mnt/lustre/groups/CBBI0818/TZ/data/qc_results/step5_final/qc_final.fam", header=FALSE)
colnames(fam) <- c("FID", "IID", "Father", "Mother", "Sex", "Phenotype")

# Merger
covariates <- merge(pca[, c("FID", "IID", paste0("PC", 1:10))], 
                     fam[, c("FID", "IID", "Sex")], 
                     by=c("FID", "IID"))

# Sauvegarder (format: FID IID Sex PC1 PC2 ... PC10)
write.table(covariates, "/mnt/lustre/groups/CBBI0818/TZ/data/qc_results/step5_final/covariates.txt",
            quote=FALSE, row.names=FALSE, col.names=TRUE)

cat("Fichier de covariables créé avec succès\n")
cat("Colonnes: FID IID Sex PC1-PC10\n")
RSCRIPT

chmod +x ${OUTPUT_DIR}/step5_final/create_covariates.R
Rscript ${OUTPUT_DIR}/step5_final/create_covariates.R

################################################################################
# ÉTAPE 6: RAPPORT FINAL
################################################################################
echo ""
echo "[ÉTAPE 6] Génération du rapport final..."

cat > ${OUTPUT_DIR}/QC_REPORT.txt << REPORT
================================================================================
                    RAPPORT DE CONTRÔLE QUALITÉ GWAS
================================================================================
Date de l'analyse: $(date)
Pipeline: QC Robuste pour GWAS + PRS

--------------------------------------------------------------------------------
STATISTIQUES INITIALES
--------------------------------------------------------------------------------
$(awk 'NR==2 {print "Nombre de SNPs: " $4}' ${OUTPUT_DIR}/step1_initial/initial_stats.log)
$(awk 'NR==2 {print "Nombre d'individus: " $2}' ${OUTPUT_DIR}/step1_initial/initial_stats.log)

--------------------------------------------------------------------------------
STATISTIQUES APRÈS QC
--------------------------------------------------------------------------------
$(awk 'END {print "SNPs après QC: " $4}' ${OUTPUT_DIR}/step5_final/final_stats.log)
$(awk 'END {print "Individus après QC: " $2}' ${OUTPUT_DIR}/step5_final/final_stats.log)

--------------------------------------------------------------------------------
FILTRES APPLIQUÉS
--------------------------------------------------------------------------------
SNPs:
  - Call rate: > $(echo "1 - ${GENO_SNP}" | bc)
  - MAF: > ${MAF_THRESHOLD}
  - HWE p-value: > ${HWE_PVALUE}

Échantillons:
  - Call rate: > $(echo "1 - ${GENO_INDIV}" | bc)
  - Discordances de sexe: retiré
  - Outliers hétérozygotie (±${HET_SD} SD): retiré
  - Individuals apparentés (PI_HAT > ${IBD_THRESHOLD}): retiré
  - Outliers de population (±6 SD sur PCs): retiré

--------------------------------------------------------------------------------
FICHIERS DE SORTIE PRINCIPAUX
--------------------------------------------------------------------------------
Données génotypiques finales:
  - ${OUTPUT_DIR}/step5_final/${FINAL_PREFIX}.bed/bim/fam

Covariables pour GWAS:
  - ${OUTPUT_DIR}/step5_final/covariates.txt
    (contient: FID, IID, Sex, PC1-PC10)

Composantes principales:
  - ${OUTPUT_DIR}/step4_population/pca.eigenvec (20 PCs)
  - ${OUTPUT_DIR}/step4_population/pca.eigenval

Graphiques QC:
  - ${OUTPUT_DIR}/plots/heterozygosity_distribution.png
  - ${OUTPUT_DIR}/plots/pca_pc1_pc2.png
  - ${OUTPUT_DIR}/plots/pca_pc1_pc3.png
  - ${OUTPUT_DIR}/plots/pca_scree.png

--------------------------------------------------------------------------------
PROCHAINES ÉTAPES RECOMMANDÉES
--------------------------------------------------------------------------------
1. GWAS avec EMMAX:
   - Utiliser: ${OUTPUT_DIR}/step5_final/${FINAL_PREFIX}.bed/bim/fam
   - Covariables: ${OUTPUT_DIR}/step5_final/covariates.txt
   - Ajuster pour sexe + 10 PCs

2. Pour PRS de qualité:
   - Vérifier la concordance d'ancestralité avec cohorte de base PRS
   - Utiliser les mêmes SNPs pour le PRS
   - Ajuster pour les mêmes covariables (sexe + PCs)

3. Vérifications additionnelles recommandées:
   - Vérifier inflation génomique (lambda) après GWAS
   - QQ-plot des p-values GWAS
   - Manhattan plot pour identifier signaux

================================================================================
REPORT

echo ""
echo "=================================="
echo "Pipeline QC terminé avec succès!"
echo "=================================="
echo ""
echo "Consultez le rapport complet: ${OUTPUT_DIR}/QC_REPORT.txt"
echo ""
echo "Fichiers principaux:"
echo "  - Données: ${OUTPUT_DIR}/step5_final/${FINAL_PREFIX}.bed/bim/fam"
echo "  - Covariables: ${OUTPUT_DIR}/step5_final/covariates.txt"
echo "  - PCA: ${OUTPUT_DIR}/step4_population/pca.eigenvec"
echo ""
echo "Nombre final d'échantillons:"
wc -l < ${OUTPUT_DIR}/step5_final/${FINAL_PREFIX}.fam
echo ""
echo "Nombre final de SNPs:"
wc -l < ${OUTPUT_DIR}/step5_final/${FINAL_PREFIX}.bim
echo ""

################################################################################
# FIN DU PIPELINE
################################################################################
