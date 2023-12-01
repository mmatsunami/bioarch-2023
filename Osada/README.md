# クオリティチェック，マッピング，バリアントコール，フィルタリング

## クオリティチェック

```
alias fastp='singularity exec /usr/local/biotools/f/fastp\:0.23.4--hadf994f_2 fastp'
alias samtools='singularity exec /usr/local/biotools/s/samtools\:0.1.19--h94a8ba4_6 samtools'
alias bwa='singularity exec /usr/local/biotools/b/bwa\:0.7.8--hed695b0_5 bwa'
alias samblaster='singularity exec /usr/local/biotools/s/samblaster\:0.1.26--hc9558a2_0 samblaster'
alias gatk='singularity exec /usr/local/biotools/g/gatk4\:4.4.0.0--py36hdfd78af_0 gatk'
alias bcftools='singularity exec /usr/local/biotools/b/bcftools\:1.18--h8b25389_0 bcftools'
alias bedtools='singularity exec /usr/local/biotools/b/bedtools\:2.31.1--hf5e1c6e_0 bedtools'
```
### fastqファイルの確認
次の2つのfastqファイルを使います．R1が名前についているファイルはforwardリード，R2が付いているファイルはreverseリードです．長さはそれぞれ150bpです．
```
/home/bioarchaeology-pg/data/3/SP01_R1.fq.gz
/home/bioarchaeology-pg/data/3/SP01_R2.fq.gz
```
まずはファイルを確認してみましょう．ファイルはgzipで圧縮されていますが，新しめのシェルであればlessコマンドで直接中身を見ることができます．
```
less /home/bioarchaeology-pg/data/3/SP01_R1.fq.gz
```
中身を開いてみると，4行で1つのリードを表したデータを見ることができます．4行目はクオリティスコアです．
