import soundfile as sf
import neural_inference.ni as ni

wavfile1 = "sample_audios/Andrew1.wav"
wavfile2 = "sample_audios/Andrew2.wav"


y1, sr1 = ni.read(wavfile1)
y2, sr2 = ni.read(wavfile2)

pp1=ni.embeddings(y1, sr1)
pp2=ni.embeddings(y2, sr2)

dtw = ni.dtw(pp1,pp2)