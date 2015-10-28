# call_marking_gui
GUI tool for marking calls, designed with heavy manual operation

This is a Matlab GUI tool for marking calls in recordings from animals. The program can take multi-channel recordings and allow user to switch between different channels. When a call is marked the channel in which the marking occurs is saved.

Input to the program should be a `MAT` file containing `sig` (sound recording) and `fs` (sampling frequency in Hz). Below is an example of a file with 32 channel, 4 seconds of recording, sampled at 250 kHz:
```
Name         Size          Bytes       Class     Attributes
  fs           1x1              8      double              
  sig    1000002x32     256000512      double        
```

Operation:
* **Load file**  Load data file in `MAT` format
* **Detect calls**  Initital simple energy detection of calls and display the spectrogram and time seris
* **Add mark**  Add a call mark on where you click
* **Remove mark**  Remove all marks enclosed in a dragged rectangle
* **Done & Save**  Save the marked calls in a new file with default filename: `xxx_detect.mat`, where `xxx` is the filename of the loaded data file. Therefore, the call marking results and the raw data are kept separate.
* Note the `MATLAB` figure tools `Zoom in`, `Zoom out`, and `Pan` are active in this program. Scroll bar on the right-hand side of the spectrogram is used to control the coloraxis of spectrogram.
* After saving `xxx_detect.mat`, you can use **Load file** to load either `xxx_detect.mat` or `xxx.mat` to review what has been saved.

