# PCR_recipe
PCRで必要な各種試薬の量やTm値から求めるthermal cyclerの温度設定を自動的に出力するソフト。(required python > 3.7)

## 必要なもの
+ dataclasses
+ Biopython
+ pandas
+ numpy
いずれも2020-03-29現在で最新版を使用している。pip or anacondaで適宜入れて下さい。

## 使用方法
1. PCR GUI.pyとpcr_recipe.pyを同じディレクトリに格納する。
2. PCR GUI.pyを走らせる
3. 以下の画面が出てくるので、primer、プライマー濃度、増幅領域、テンプレート濃度 etc. を入力
   (使うメーカーは、2020-03-29現在で**KOD or PrimeSTARのみ**対応している。KOD FXとか細かいものは適宜増やしていく予定。)
![image](https://user-images.githubusercontent.com/41857834/112939815-167ecd80-9167-11eb-98c5-5e1a57ecbd66.png)
4. 次へをクリックすれば、Thermal cyclerの設定と必要試薬量が一覧で見られる。また、「HTMLに出力」をクリックすると1. のディレクトリに
HTMLファイルが保存されるので、ブラウザで印刷することも可能。
![image](https://user-images.githubusercontent.com/41857834/112939853-29919d80-9167-11eb-94f7-eebc2def2974.png)
