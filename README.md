# Spectrogram reassignment

This function computes reassigned version of the spectrogram. The algorithm is based on [1], some parts are adopted from [2]. The idea is to first compute conventional spectrogram, then find optimal (in a sense of energy) time and frequency positions and reassigns values in the spectrogram to this new positions.

You can provide additional properties to the function. Detailed description of the properties are below and the full list with possible values is in [Syntax section](##Syntax).

You have a choice between reassigning spectrogram either as short-time Fourier transform (STFT) or power spectrum density (STFT normalized by window energy and sampling rate). Controlled by 'psd' property.

You can pre-process signal by padding it. Several types of padding are available: with zeros, with a constant, periodical and symmetrical padding. This can be done by providing a type of padding in 'pad' property.

Sometimes, extreme values (outliers) appear on a spectrogram (reassigned and original), for instance, when using too many overlapping points between the segments. You can 'crop' them by providing to a desired value of percentile of the data, above which all the values will be set to the highest values. It works similarly to thresholding. To use this feature specify percentile as 'crop' property.

It is possible to obtain spectrogram on finer grid with this method. For that, you have to specify in the properties three additional parameters: new steps or new dimensions of time and frequency, specify if it is new steps or dimensions. Usage of this feature is not recommended, because even though it is possible to increase number of points in the output matrix the original number of points is not changing. So, you might obtain very sparse matrix. To partially solve it you can use interpolation (set 'interp' parameter), which interpolates additional points with values of its neighbouring points. A use of this feature is also not advised, because it decreases the sharpness of reassigned spectrogram and contradicts the whole idea of the reassignment.


# Syntax
REASSPECGRAM spectrogram reassignment
 RS = REASSPECGRAM(SIG,WIN,OVLAP,NFFT,FS,'PROPERTYNAME',PROPERTYVALUE)
 RS = REASSPECGRAM(SIG,WIN,OVLAP,NFFT,FS,PROPERTYSTRUCTURE)
[RS,S] = REASSPECGRAM(...)
[RS,FNEW,TNEW] = REASSPECGRAM(...)
[RS,FNEW,TNEW,S,FORIG,TORIG] = REASSPECGRAM(...)

|Input  |                   |
|-------|-------------------|
|SIG    | analysed signal of length N points. |
|WIN    | smoothing window. If WIN is a vector, the signal SIG is divided into segments of length equal to the length of the window, and each segment is filtered with WIN. If WIN is an integer, then SIG is divided into segments of length equal to the specified integer and filtered with Hamming window. If no window is provided, then Hamming window of length N/10 is used by default. |
|OVLAP  | number of overlapping points between two adjacent windows. Is used to compensate energy loss at the ends of smoothing windows. Should be specified as an integer smaller than the length of the window. If no OVLAP is given, the default value is used to obtain a 50% overlap. |
|NFFT   | number of Fourier transform points. If no input is provided, the default value is chosen as next power of 2 greater then the length of the window. |
|FS     | sampling frequency in Hz.|

|Output |                   |
|-------|-------------------|
| RS    | A matrix with reassigned spectrogram. Here, rows are frequencies and columns are time points.|
| S     | A matrix with original (not reassigned) spectrogram. Here, rows are frequencies and columns are time points. |
| FNEW  | Frequency vector corresponding to rows in reassigned spectrogram matrix. |
| TNEW  | Time vector corresponding to columns in reassigned spectrogram matrix. |
| FORIG | Frequency vector corresponding to rows in spectrogram matrix. |
| TORIG | Time vector corresponding to columns in spectrogram matrix. |

|Property  | Possible values      |                   |
|----------|----------------------|-------------------|
| 'psd'    | true, **false** | If true, power spectral density (PSD) is reassigned. If not, squared absolute value of short-time Fourier transform is used. |
| 'pad'    | **'zeros'** , 'const' , 'periodic' , symmetric' | If specified before processing a signal will be padded according to the chosen method. Padding can be done with zeros ('zeros'); with constant  equal to the first and last values of the signal ('const'); 'periodic' continues signal in both directions periodically; 'symmetric' reflects signal symmetrically on both sides. If padding is chosen it is added on both slides, the length of the padding is half of the window length. |
| 'crop'   | true, **false** | Sometimes there might be extreme values in the final time-frequency representation. You can "crop" such values, by specifying percentile, above which the values should not be considered (analogy with outliers), value of 99.5-99.9% is reasonable. |
| 'size'   | 2-element vector| You can change the size of the output TFR matrix.  Specify a new size as two-element vector: number of rows (frequency) and number of columns (time). Not recommended, because the number of elements in the TFR matrix is fixed. |
| 'step'   | 2-element vector| Similar to 'size' parameter, but instead of specifying new size you provide new frequency and time sampling. |
| 'interp' | true, **false** | Interpolate neighbouring points. Can be useful when using new size/sampling, because it increases the number of points.  Otherwise, it is not recommended to use it, because the idea of reassignment is lost like that. |
By default all properties are unset. Property can be passed as name-value pair
or as a structure. If field requires logical true/false input, every non-empty,
finite and not 'NaN' input is treated as logical true, e.g. true, 'yes', 1
etc. will be validated as logical true.

