# クオリティチェック，マッピング，バリアントコール，フィルタリング

## 準備
データ解析をするためのワーキングディレクトリを作っておきます．以後，このディレクトリ下で作業を行いましょう．
```bash
cd #自分のホームディレクトリに移動する
mkdir working
cd working
```
`mkdir`コマンドでフォルダを作り，`cd`コマンドでその場所で移動します．
## 実習で用いるデータ
実習で用いるデータは，日本列島で採取された架空の2倍体生物，_Fictus yaponesiae_ から得られたショートリード配列です．日本の5都市（札幌，SP；仙台，SD；東京，TK；大阪，OS；福岡，FK）から10個体ずつサンプリングされています．ゲノム配列は小さく，約23Mbです．
## 使用するプログラム
以下のコマンドをコピーしてターミナル上にペーストすることによって，インストール無しでプログラムを使うことができます．
```bash
alias fastp='singularity exec -B /lustre8/home,/home /usr/local/biotools/f/fastp\:0.23.4--hadf994f_2 fastp'
alias samtools='singularity exec -B /lustre8/home,/home /usr/local/biotools/s/samtools:1.18--hd87286a_0 samtools'
alias bwa='singularity exec -B /lustre8/home,/home /usr/local/biotools/b/bwa\:0.7.8--hed695b0_5 bwa'
alias samblaster='singularity exec -B /lustre8/home,/home /usr/local/biotools/s/samblaster\:0.1.26--hc9558a2_0 samblaster'
alias gatk='singularity exec -B /lustre8/home,/home /usr/local/biotools/g/gatk4\:4.4.0.0--py36hdfd78af_0 gatk'
alias bcftools='singularity exec -B /lustre8/home,/home /usr/local/biotools/b/bcftools\:1.18--h8b25389_0 bcftools'
alias bedtools='singularity exec -B /lustre8/home,/home /usr/local/biotools/b/bedtools\:2.31.1--hf5e1c6e_0 bedtools'
alias picard='singularity exec -B /lustre8/home,/home /usr/local/biotools/p/picard\:3.0.0--hdfd78af_1 picard'
```
試しに次のコマンドを売ってみましょう．オプションの一覧が表示されれば成功です．
```bash
samtools
```
## クオリティチェック（quality checking）
### fastqファイルの確認
次の2つのfastqファイルを使います．R1が名前についているファイルはforwardリード，R2が付いているファイルはreverseリードです．長さはそれぞれ150bpです．
```
/home/bioarchaeology-pg/data/Osada/SP01_R1.fq.gz
/home/bioarchaeology-pg/data/Osada/SP01_R2.fq.gz
```
これらのファイルを参照したいのですが，bioarhaeology-pgのような文字列をいちいち打つのも面倒なので，シンボリックリンクを作ってみます．
```bash
ln -s /home/bioarchaeology-pg/data/Osada Osada
```
このようにすることで，working/Osadaが/home/bioarchaeology-pg/data/Osadaを指すようになります．

まずはファイルを確認してみましょう．ファイルはgzipで圧縮されていますが，新しめのシェルであれば`less`コマンドで直接中身を見ることができます．
```bash
less Osada/SP01_R1.fastq.gz
```
中身を開いてみると，4行で1つのリードを表したデータを見ることができます．4行目はクオリティスコアです．`less`を終了するには`q`を押します．

ファイルの行数を数えてみましょう．`zcat`コマンドでgzip圧縮されたファイルを標準出力（スクリーン）に表示することができます．その結果を`wc`コマンドにパイプして行数を数えることができます．パイプは`|`という文字で表され，標準出力を次のコマンドの引数として渡すためのものです．
```bash
zcat Osada/SP01_R1.fastq.gz | wc
```
次のような出力があるはずです．
```
10345564 10345564 831706850
```
左から，行数，単語数，バイト数です．fastq配列にはスペースが入っていませんから行数=単語数になります．

行数を4で割れば配列の本数になります．別の方法として，`grep`コマンドでファイル内の`@`の数を数えることも可能です．`grep`は入力文字列の検索を行うことができ，非常によく使うコマンドです．
```bash
zcat Osada/SP01_R1.fastq.gz | grep '@' | wc
```
同様に，reverseリードが格納されている`SP01_R2.fq.gz`ファイルについても行数を表示します．bwaなどのマッピングソフトは，fastqファイルが2つに分かれている場合には，上から1つずつ配列を読み込んでいきそれらがペアであるという前提で動いていきます．したがって，`_R1.fastq.gz`と`_R2.fastq.gz`には同じ数の配列が入っている必要があります．2つのファイルで配列の本数が違う場合はSeqkitなどのソフトウェアを用いて配列をきれいに整理してあげる必要があります．

fastpソフトウェアを使ってfastqファイルのフィルタリングをしてみます．標準的なアダプター配列の除去も行ってくれます．次のコマンドは`fastp`プログラムを用いてペアエンド配列を解析し，フィルターした配列を別のファイル`SP01_R1.clean.fastq.gz`と`SP01_R1.clean.fastq.gz`に出力します．
```bash
fastp -i Osada/SP01_R1.fastq.gz -I Osada/SP01_R2.fastq.gz -o SP01_R1.clean.fastq.gz -O SP01_R2.clean.fastq.gz
```
ワーキングディレクトリに`fastp.html`というファイルが出力されるので確認してみましょう．sftpソフトを用いて自分のPCにファイルをダウンロードし，ブラウザで確認することができます．
## マッピング（mapping）
### fastqファイルの縮小
fastqファイルの確認とフィルタリングが終わったら，マッピングを行います．今回は練習ですので，解析するリード数を大幅に減らしてみましょう．
```bash
zcat SP01_R1.clean.fastq.gz | head -1000000 | gzip > SP01_R1.reduced.fastq.gz
zcat SP01_R2.clean.fastq.gz | head -1000000 | gzip > SP01_R2.reduced.fastq.gz
```
`head`コマンドは引数の行数だけ頭から取り出して標準出力に送ります．この例ではその出力をパイプで`gzip`に渡してgzファイルに圧縮しています．それぞれのfastqファイルは25万配列が入っています．
### bwaによるマッピング
今回は，最もよく使われているマッピングソフトウェア，bwaを利用します．まずは次のコマンドでワーキングディレクトリにゲノムのリファレンス配列をコピーします．
```bash
cp Osada/yaponesia_genome.fasta ./
```
bwaでインデックスファイルを作成します．この作業により，相同性検索の高速化が可能です．
```bash
bwa index yaponesia_genome.fasta
```
マッピングを行うと，samフォーマットのアラインメントが標準出力に出力されます．ここでは`bwa mem`を使ってペアエンド配列をマッピングします．リファレンス配列，リード配列1，リード配列2の順番で入力します．ここでは，オプション`-R`を使って`-R "@RG\tID:SP01\tLB:LB\tSM:SP01\tPL:ILLUMINA"`のようにリードグループを指定してマッピングを行っています．これは後で使うGATKで必要になるもので，異なったサンプルやライブラリから由来するリードを区別するためのタグになります．他のソフトウェアでもこの情報が使われる場合があります．
```bash
bwa mem -R "@RG\tID:SP01\tLB:LB\tSM:SP01\tPL:ILLUMINA" yaponesia_reference.fasta SP01_R1.clean.fastq.gz SP01_R2.clean.fastq.gz
```
>[!NOTE]
>ものすごい勢いで出力が画面に流れるので`Ctrl+C`でプログラムを止める必要があります．

同時に，Samblasterソフトウェアを使ってPCR duplicatesをマークします．Samblasterはsamフォーマットの入力を受け取り，マーク済みのsamフォーマットを標準出力に返します．結果をファイルとして保存するにはsamblasterの出力をsamtoolsで受け取って，ファイルに出力します．コマンドは次のようになります．
```bash
bwa mem -R "@RG\tID:SP01\tLB:LB\tSM:SP01\tPL:ILLUMINA" yaponesia_reference.fasta SP01_R1.reduced.fastq.gz SP01_R2.reduced.fastq.gz | samblaster | samtools sort -O BAM -o SP01.temp.bam
```
一番最後のプロセスをよく見てみましょう．samtoolsはsamフォーマットのファイルを扱うためのスタンダードなソフトウェアです．受け取ったsamフォーマットのファイルを`samtools sort`コマンドでリファレンスゲノム配列の座標の順番に並べ直します．この作業はほとんどのバリアントコールの前に必要な作業です．ただし，ここでは`-o`オプションにより出力ファイルが指定されているので，標準出力ではなくファイルに出力されます．`-O BAM`はbamフォーマットで出力することを指定しています．大規模なプロジェクトであればファイルサイズを減らすために`-O CRAM`でcramファイルとして出力することも考えてみましょう．

出来上がったbamファイルは次のコマンドで確認できます．samファイルは`less`コマンドで直接表示できますが，bamファイルは圧縮されたファイルなので`samtools view`を使って表示する必要があります．表示が流れるのを防ぐために`less`コマンドにパイプします．
```bash
samtools view SP01.temp.bam | less
```
cramファイルも同様に表示できますが，cramファイルを扱うときには常にオプション`-t`でマッピングに使ったリファレンス配列を指定する必要があります．

今後の作業のためにbamファイルのインデックスを作っておく必要があります．
```bash
samtools index SP01.temp.bam
```
indexを作っておくと，次のようなコマンドを用いてbamファイルの一部を高速に表示することが可能です．
```bash
#　1番染色体の1000～2000番目の塩基にマッピングされるリードを表示
samtools view SP01.temp.bam chr1:1000-2000
```
マッピングの結果をPicardソフトウェアで調べてみましょう．`picard CollectWgsMetircs`で平均カバー率などを見ることができます．また，複数サンプルの結果をMultiQCソフトウェアでまとめて表示することも可能です．
```bash
picard CollectWgsMetrics -I SP01.bam -O SP01.picardCWM.txt -R yaponesia_reference.fasta
```
結果は指定した`SP01.picardCWM.txt`に記録されます．結果を`less`コマンドで見てみましょう．MEAN_COVERAGEやMEDIAN＿COVERAGE，カバー率のヒストグラムデータを見ることができます．
## バリアントコール（variant calling）
本実習ではバリアントコールにはGATKを使います．VCFフォーマットを出力するバリアントコール，GVCFフォーマットを出力するバリアントコールの両方を試してみましょう．どちらの場合も`gatk HaplotypeCaller`を使用します．

HaoplotypeCallerを行うには，picardを用いてリファレンス配列のdictファイルとfaiファイルを作る必要があります．どちらもインデックスファイルです．
```bash
picard CreateSequenceDictionary -R yaponesia_reference.fasta -O yaponesia_reference.dict
samtools faidx yaponesia_reference.fasta
```


まずはVCFファイルを出力します．bamファイルは先ほど作成した縮小版ではなく，すべてのリードをマッピングしたものを使います．時間を節約するために解析する領域を`-L`オプションで限定します．
```bash
gatk HaplotypeCaller -I Osada/SP01.bam -R yaponesia_reference.fasta -O SP01.vcf.gz -L chr1:1-1000000
```
結果を画面で確認します．
```bash
less SP01.vcf.gz
```
次に，GVCFフォーマットの出力を行います．
```bash
gatk HaplotypeCaller -I SP01.bam -R yaponesia_reference.fasta -O SP01.gvcf.gz -L chr1:1-1000000 -ERC GVCF
```
結果を画面で確認します．
```bash
less SP01.gvcf.gz
```
GVCFフォーマットには，変異のない部分の情報が含まれていることがわかります．
>[!NOTE]
>現在のバージョンのGATKでは，カバー率（DP）が0のところも`0/0`とコールされる仕様となっていますので注意が必要です．
### 多検体の処理
>[!NOTE]
>多数のサンプルのバリアントコールを行う場合は，シェルスクリプトを作成すると効率よく進みます．より簡便な方法では`xargs`や`paralllel`コマンドも便利です．
>
>以下は多サンプルを効率よくマッピングするシェルスクリプトの例です．本実習では時間の都合上実際に動かすことはありませんが，興味のある方はフォローしてください．
>
>まず，サンプル名が1行に1サンプルずつ書かれたテキストファイルを用意します．次は`samples.txt`の中身です．
>```
># samples.txt
>SP01
>TK01
>```
>次のようなシェルスクリプトを作り，`mapping.sh`として保存します（ファイルはディレクトリ`Osada`にもあります）．
>```bash
># mapping_calling.sh
>#! /bin/bash
>
># ここはNIGスパコンで使うためのおまじないなので気にしない
>alias samtools='singularity exec -B /lustre8/home,/home /usr/local/biotools/s/samtools:1.18--hd87286a_0 samtools'
>alias bwa='singularity exec -B /lustre8/home,/home /usr/local/biotools/b/bwa\:0.7.8--hed695b0_5 bwa'
>alias samblaster='singularity exec -B /lustre8/home,/home /usr/local/biotools/s/samblaster\:0.1.26--hc9558a2_0 samblaster'
>alias gatk='singularity exec -B /lustre8/home,/home /usr/local/biotools/g/gatk4\:4.4.0.0--py36hdfd78af_0 gatk'
>
># ここからが本体
>while read LINE; do
>  ID=$LINE
>  echo $ID
>  RG="@RG\tID:${ID}\tLB:LB\tSM:${ID}\tPL:ILLUMINA"
>  bwa mem -R $RG yaponesia_reference.fasta ${ID}_R1.fastq.gz ${ID}_R2.fastq.gz | samblaster | samtools sort -O BAM -o ${ID}.bam
>  samtools index ${ID}.bam
>  gatk HaplotypeCaller -I ${ID}.bam -R yaponesia_reference.fasta -O ${ID}.gvcf.gz -ERC GVCF 
>done < samples.txt
>```
>この例では，`sample.txt`から受け取った名前を変数`$ID`に格納し，マッピング，重複リードの除去，ソート，バリアントコールによるGVCFファイルの生成までをひとまとめに行っています．実際の解析では，bamファイルやgvcfファイルはそれぞれ異なったディレクトリに保存しておくと良いでしょう．
>
>シェルスクリプト`mapping_calling.sh`をワーキングディレクトリにコピーしてから実行することによって，2つのGVCFファイル，`SP01.gvcf.gz`, `TK01.gvcf.gz`がワーキングディレクトリに作成されます．`Osada`ディレクトリにも同じものが入っています．
>2つのGVCFファイルをGATKを用いて合体させるには，以下のコマンドを用います．
>```bash
>gatk GenomicsDBImport -V SP01.gvcf.gz -V TK01.gvcf.gz --genomicsdb-workspace-path yaponesiadb -L chr1
>gatk GenotypeGVCFs -R yaponesia_reference.fasta -V gendb://yaponesiadb -O SP01_TK01.vcf.gz
>```
>たくさんのサンプルがある場合は，サンプルマップと呼ばれる，サンプル名とGVCFファイル名を対応付けたテキストファイルを作成すると便利でしょう．[GATKウェブサイト](https://gatk.broadinstitute.org/hc/en-us/articles/360036883491-GenomicsDBImport)
>一度データベースとしてGVCFファイルの情報をひとまとめにしてから，`GATK GenotypeGVCFs`というコマンドでジョイントコールを行います．この練習データにはindelの変異は入っていませんので，検出される変異はすべてSNVになります．
## フィルタリング（filtering）
この実習では簡単にハードフィルタをかけてみます．GATKでSNVのために推奨されているハードフィルターを描けます．
```bash
gatk VariantFiltration \
  -V Osada/SP01_TK01.vcf.gz \
  -filter "QD<2.0" --filter-name "QD2" \
  -filter "QUAL<30.0" --filter-name "QUAL30" \
  -filter "SOR>3.0" --filter-name "SOR3" \
  -filter "FS>60.0" --filter-name "FS60" \
  -filter "MQ<40.0" --filter-name "MQ40" \
  -filter "MQRankSum<-12.5" --filter-name "MQRankSum-12.5" \
  -filter "ReadPosRankSum<-8.0" --filter-name "ReadPosRankSum-8" \
  -filter "ExcessHet>30.0" --filter-name "ExcessHet30" \
  -O SP01_TK01.Hardfilter.vcf.gz \
  -verbosity ERROR
```
それぞれのSNVについて，`INFO`列をチェックし，`-filter`で指定されたフィルターに引っかかるものがあると，それに対応するフィルター名`--filter-name`を`FILTER`列に書き出します．この時点ではSNVは取り除かれていません．それぞれのフィルターの意味については[GATKウェブサイト](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants)を参照してください．また，ここでは`ExcessHet`として示されているヘテロ接合過剰に関するフィルターもかけています．

結果を表示するといくつかの変異がフィルターに引っかかったのがわかります．フィルターをパスした変異には`PASS`という文字が`FILTER`列に記入されているのがわかります．

最後に，フィルターに引っかかった変異を取り除いてvcfファイルを完成させます．
```bash
gatk SelectVariants -V SP01_TK01.Hardfilter.vcf.gz -O SP01_TK01.Hardfiltered.vcf.gz --exclude-filtered
```
以下のコマンドでSNVの数を簡便にカウントできます．`grep -v "^#"`は先頭に`#`を含むヘッダー行以外を表示する方法です．
```bash
# フィルター前の変異数をカウント
zcat SP01_TK01.Hardfilter.vcf.gz | grep -v "^#" | wc
# フィルター後の変異数をカウント
zcat SP01_TK01.Hardfiltered.vcf.gz | grep -v "^#" | wc
```
いくつかのSNVがフィルタリングされたことがわかります．

vcfファイルの操作には`bcftools`が便利です．次のコマンドは，`SP01_TK01.Hardfiltered.vcf.gz`から，サンプルTK01だけを抜き出します．
```bash
bcftools view -s TK01 -Oz -o SP01.Hardfiltered.vcf.gz SP01_TK01.Hardfiltered.vcf.gz
```
他にも便利な機能がありますので，[bcftoolsのマニュアル](https://samtools.github.io/bcftools/bcftools.html)を確認すると良いでしょう．
