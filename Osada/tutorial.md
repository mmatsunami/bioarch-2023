# クオリティチェック，マッピング，バリアントコール，フィルタリング

## 準備
データ解析をするためのワーキングディレクトリを作っておきます．以後，このディレクトリ下で作業を行いましょう．
```
%mkdir working
%cd working
```
`mkdir`コマンドでフォルダを作り，`cd`コマンドでその場所で移動します．`%`はプロンプトを表しているので，入力する必要はありません．ほかにも`>`や`$`がプロンプトを表す文字として，このようなテキストでは用いられています．
## 実習で用いるデータ
実習で用いるデータは，日本列島で採取された架空の生物，_Fictus yaponesiae_ から得られたショートリード配列です．日本の5都市（札幌，SP；仙台，SD；東京，TK；大阪，OS；福岡，FK）から10個体ずつサンプリングされています．ゲノム配列は小さく，約100Mbです．
## クオリティチェック（quality checking）

```
alias fastp='singularity exec -B /lustre8/home,/home /usr/local/biotools/f/fastp\:0.23.4--hadf994f_2 fastp'
alias samtools='singularity exec -B /lustre8/home,/home /usr/local/biotools/s/samtools\:0.1.19--h94a8ba4_6 samtools'
alias bwa='singularity exec -B /lustre8/home,/home /usr/local/biotools/b/bwa\:0.7.8--hed695b0_5 bwa'
alias samblaster='singularity exec -B /lustre8/home,/home /usr/local/biotools/s/samblaster\:0.1.26--hc9558a2_0 samblaster'
alias gatk='singularity exec -B /lustre8/home,/home /usr/local/biotools/g/gatk4\:4.4.0.0--py36hdfd78af_0 gatk'
alias bcftools='singularity exec -B /lustre8/home,/home /usr/local/biotools/b/bcftools\:1.18--h8b25389_0 bcftools'
alias bedtools='singularity exec -B /lustre8/home,/home /usr/local/biotools/b/bedtools\:2.31.1--hf5e1c6e_0 bedtools'
```
### fastqファイルの確認
次の2つのfastqファイルを使います．R1が名前についているファイルはforwardリード，R2が付いているファイルはreverseリードです．長さはそれぞれ150bpです．
```
/home/bioarchaeology-pg/data/Osada/SP01_R1.fq.gz
/home/bioarchaeology-pg/data/Osada/SP01_R2.fq.gz
```
これらのファイルを参照したいのですが，bioarhaeology-pgのような文字列をいちいち打つのも面倒なので，シンボリックリンクを作ってみます．
```
%ln -s /home/bioarchaeology-pg/data/Osada Osada
```
このようにすることで，working/Osadaが/home/bioarchaeology-pg/data/Osadaを指すようになります．

まずはファイルを確認してみましょう．ファイルはgzipで圧縮されていますが，新しめのシェルであれば`less`コマンドで直接中身を見ることができます．
```
%less SP01_R1.fq.gz
```
中身を開いてみると，4行で1つのリードを表したデータを見ることができます．4行目はクオリティスコアです．

ファイルの行数を数えてみましょう．`zcat`コマンドでgzip圧縮されたファイルを標準出力（スクリーン）に表示することができます．その結果を`wc`コマンドにパイプして行数を数えることができます．パイプは`|`という文字で表され，標準出力を次のコマンドの引数として渡すためのものです．
```
%zcat SP01_R1.fq.gz | wc
```
行数を4で割れば配列の本数になります．別の方法として，`grep`コマンドでファイル内の`@`の数を数えることも可能です．`grep`は入力文字列の検索を行うことができ，非常によく使うコマンドです．
```
%zcat SP01_R1.fastq.gz | grep '@' SP01_R1.fq.gz | wc
```
同様に，reverseリードが格納されている`SP01_R2.fq.gz`ファイルについても行数を表示します．bwaなどのマッピングソフトは，fastqファイルが2つに分かれている場合には，上から1つずつ配列を読み込んでいきそれらがペアであるという前提で動いていきます．したがって，`_R1.fastq.gz`と`_R2.fastq.gz`には同じ数の配列が入っている必要があります．2つのファイルで配列の本数が違う場合はSeqkitなどのソフトウェアを用いて配列をきれいに整理してあげる必要があります．

fastpソフトウェアを使ってfastqファイルのフィルタリングをしてみます．標準的なアダプター配列の除去も行ってくれます．次のコマンドは`fastp`プログラムを用いてペアエンド配列を解析し，フィルターした配列を別のファイル`SP01_R1.clean.fastq.gz`と`SP01_R1.clean.fastq.gz`に出力します．
```
%fastp -i SP01_R1.fastq.gz -I SP01_R1.fastq.gz -o SP01_R1.clean.fastq.gz -O SP01_R2.clean.fastq.gz
```
ワーキングディレクトリに`fastp.html`というファイルが出力されるので確認してみよう．sftpソフトを用いて自分のPCにファイルをダウンロードし，ブラウザで確認することができる．
### マッピング
fastqファイルの確認とフィルタリングが終わったら，マッピングを行う．今回は，最もよく使われているマッピングソフトウェア，bwaを利用する．

まずは次のコマンドでワーキングディレクトリにゲノムのリファレンス配列をコピーする．
```
%cp /home/bioarchaeology-pg/data/Osada/yaponesia_genome.fasta ./
```
bwaでインデックスファイルを作成します．この作業により，相同性検索の高速化が可能です．
```
%bwa index yaponesia_genome.fasta
```
マッピングを行うと，samフォーマットのアラインメントが標準出力に出力されます．ここでは`bwa mem`を使ってペアエンド配列をマッピングします．リファレンス配列，リード配列1，リード配列2の順番で入力します．
```
%bwa mem yaponesia_reference.fasta SP01_R1.clean.fastq.gz SP01_R2.clean.fastq.gz
```




