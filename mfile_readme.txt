＜mファイル機能＞
	•	makingMat.m 実験データ(.xlsx)fileをmatファイルに切り替える
	•	angular_velocity.m 角度データより角速度データを算出する（二点近似）
	•	muscle_strength_estimation.m 生体長の時、筋の最大張力を計算する
	•	musculoskeletal6.m particle filterによる筋骨格パラメータの推定し, その結果を基にグラフ（卒論Figs.5.1-5.14）を生成
	•	musculoskeletal61.m particle filterによる筋骨格パラメータの推定（改善案）し, その結果を基にグラフ（卒論Figs.5.15-5.28）を生成
	•	result_graph.m グラフ(卒論Figs.4.8-4.17)のlayoutの修正
	•	sourceProcess_EMG.m EMGだけの信号処理（生信号→整流化→lowpassfilter）を行い、行なったデータを基にグラフ（卒論Figs.4.4-4.6）を生成
	•	sourceProcess.m 実験データ全ての信号処理（生信号→整流化→lowpassfilter）を行い、行なったデータを基にグラフ（卒論Figs.4.8-4.17）を生成

＜プログラム実行手順＞
1. makingMat.mより実験データ（experiment(19.01.30).xlsx）を（experiment.mat）に変換
2. angular_velocity.mより角速度データの算出（angular_velocity.mat）
3. sourceProcess.mより信号処理を行う（sourceProcess.mat）
4. result_graph.m ↑のデータを基にグラフ（卒論Figs.4.8-4.17）を生成
（Figs.4.8-4.12）
11行目%Process(extension); comment out
12行目   Process(flexion);
（Figs.4.13-4.17）
11行目   Process(extension); 
12行目% Process(flexion); comment out
5. sourceProcess_EMG.m EMGだけの信号処理しその結果を基にグラフ（卒論Figs.4.4-4.6）の生成
6. muscle_strength_estimation.m より筋力の推定に必要なパラメータの算出（muscle_strength_estimation.mat）
7. musculoskeletal6.m より筋骨格パラメータの推定を行う
（flexion : 45行目をコメントアウトして実行）
（extension : 44行目をコメントアウトして実行）
8. musculoskeletal61.m より改善案で筋骨格パラメータの推定を行う
（flexion : 45,106,107行目をコメントアウトして実行）
（extension : 44,104,105行目をコメントアウトして実行）

以上の順番で実行すると，
筋骨格パラメータの推定結果のmatファイルが生成（musculoskeletal6.mat）され
グラフ（卒論Figs 5.1-5.14）が生成される。
改善案で筋骨格パラメータの推定結果のmatファイルが生成（musculoskeletal61.mat）され
グラフ（卒論Figs 5.15-5.28）が生成される。

＜条件を変更する場合＞
musculoskeletal6.m, musculoskeletal61.m内の各種パラメータを変更する



