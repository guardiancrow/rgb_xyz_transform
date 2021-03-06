# 三原色と白色からRGB-XYZ変換マトリックスを求める（RGB-XYZ変換）デモ

## 概要

CIE 1931カラースペースにおいてRGBカラースペースとXYZカラースペースとの相互変換行うための変換マトリックスを求めます。

## 含まれるもの

- rgb_xyz_transform.cpp （サンプルプログラム本体）
- rgb_xyz.pdf (計算方法や補足のPDF)
- rgb_xyz.tex (PDF作成のためのtexファイル)

## 必要なもの

- boost （マトリックス変換・逆行列・極小な数 のために使用）

## 用途

ICCプロファイルや、PNG画像形式のcHRMチャンクなどで記録されている赤緑青の三原色＋白色（以下、四原色と表記）からXYZカラースペースへの変換をすることができます。またその逆もできます。

マトリックスそのものをファイルに埋め込むようなこともできますし、それを利用することもできますが、cHRMチャンクなどの場合は自分で計算する必要があります。

変換式だけを表記しているサイトは多数見かけましたが、微妙に数値が違うこと、そもそも求め方の記したサイトをほとんど見かけなかったことから公開して見ました。

## サンプルプログラムについて

サンプルプログラムではsRGBとAdobeRGB、そしてカスタム要素として赤と緑を入れ替えたsRGBのマトリックスを求めています。それぞれ違う結果であることをお確かめください。

### ビルド例

    > g++ rgb_xyz_transform.cpp -I<your_boost_path>

### サンプルで四原色のz値を指定しないのはなぜ？

```x + y + z = 1``` の関係が成り立っているのでxとyがあれば必要ないからです。

### x,y,zとX,Y,Zというように大文字小文字を書き分けてるのはなぜ？別なものなの？

はい、別のものです。小文字の方のx,y,zはCIE 1931という色空間の座標を表します。上記の通りzは省略されることが多いです。対して、大文字のX,Y,Zは三刺激値を表します。RGB色空間のR,G,Bを変換マトリックスに基づいてXYZ色空間に変換した値です。

### そもそもどうして四原色の値が違うの？

カラーマネージメントにおいて、表示するデバイス、印刷するデバイス、などごとに表示できる色が異なります。これらの差を吸収しできるだけ同じように表現しようとするために、デバイスごとに四原色を定義してあります。

また、画像やドキュメントを作成した環境と同じような色を表現するためにこれらの四原色の定義を画像ファイルなどに埋め込んである場合があります。

XYZカラースペースに変換することによってそれらの差異を「戻す」「リセットする」ことができます。そして別のデバイスに再変換することによって作成者の意図した色表現で表示することができます。

極端な話、赤色と緑色が逆に表示されるディスプレイ向けに赤色と緑色を入れ替えた四原色を指定することもできます。「そのディスプレイで作った画像を正常のディスプレイで見ても遜色ない」のようなことができてしまうのです。

現在はICCプロファイルのように、四原色の他さらにもっと詳細な情報を埋め込んである場合もありますが、基本的にCIE 1931カラースペースを介して変換することには変わりありません。

### なんでそんな面倒くさいことするの？

人間の感知できる色にはRGBカラースペースで表現できる色以外のものがあります。RGBの値がマイナスになる場合や１以上になる場合です。このような場合はXYZカラースペースで表したほうが都合がいいのです。

## その他

- 三元一次連立方程式を求めるのにガウス・ジョルダンの消去法を用いています。

## 参照

- [CIE 1931 color space : Wikipedia](https://en.wikipedia.org/wiki/CIE_1931_color_space)
- [sRGB : Wikipedia](https://en.wikipedia.org/wiki/SRGB)
- [PNG (Portable Network Graphics) Specification, Version 1.2](http://www.libpng.org/pub/png/spec/1.2/PNG-ColorAppendix.html)
- [音と色と数の散歩道](http://www.enjoy.ne.jp/~k-ichikawa/index.html)
