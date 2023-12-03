# クオリティチェック，マッピング，バリアントコール，フィルタリング

## 準備
データ解析をするためのワーキングディレクトリを作っておきます．以後，このディレクトリ下で作業を行いましょう．
```
%mkdir working
%cd working
```
`mkdir`コマンドでフォルダを作り，`cd`コマンドでその場所で移動します．%はプロンプトを表しているので，入力する必要はありません．ほかにも>や$がプロンプトを表す文字として，このようなテキストでは用いられています．
## クオリティチェック（quality checking）

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
まずはファイルを確認してみましょう．ファイルはgzipで圧縮されていますが，新しめのシェルであれば`less`コマンドで直接中身を見ることができます．
```
%less /home/bioarchaeology-pg/data/3/SP01_R1.fq.gz
```
中身を開いてみると，4行で1つのリードを表したデータを見ることができます．4行目はクオリティスコアです．

ファイルの行数を数えてみましょう．`zcat`コマンドでgzip圧縮されたファイルを標準出力（スクリーン）に表示することができます．その結果を`wc`コマンドにパイプして行数を数えることができます．パイプは`|`という文字で表され，標準出力を次のコマンドの引数として渡すためのものです．
```
%zcat /home/bioarchaeology-pg/data/3/SP01_R1.fq.gz | wc
```
行数を4で割れば配列の本数になります．別の方法として，`grep`コマンドでファイル内の`>`の数を数えることも可能です．
```
%grep '>' /home/bioarchaeology-pg/data/3/SP01_R1.fq.gz | wc
```

同様に，reverseリードが格納されている`/home/bioarchaeology-pg/data/3/SP01_R2.fq.gz`ファイルについても表示します．bwaなどのマッピングソフトは，fastqファイルが2つに分かれている場合には，上から1つずつ配列を読み込んでいきそれらがペアであるという前提で動いていきます．したがって，`_R1.fastq.gz`と`_R2.fastq.gz`には同じ数の配列が入っている必要があります．2つのファイルで配列の本数が違う場合はSeqkitなどのソフトウェアを用いて配列をきれいに整理してあげる必要があります．

fastpソフトウェアを使ってfastqファイルのフィルタリングをしてみます．標準的なアダプター配列の除去も行ってくれます．



