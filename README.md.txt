# Mini Analysis Code

The mini_analysis code can be used to analyze spontaneous events in voltage clamp recordings in Igor. This code requires "justfilters". To run mini_analysis: `mini_analysis(startwave, endwave, starttime, endtime, level, plotflag)`. This function assumes that your waves are named "sweepX" where X is the wave number.

*  `startwave` - First wave to analyze.  
*  `endwave` - Last wave to analyze.  
*  `starttime` - Time, in seconds, in each wave to start analyzing.  
*  `endtime` - Time, in seconds, in each wave to stop analyzing.  
*  `level` - Cuttoff for the slope (first derivative) to use for event detection.  
*  `plotflag` - If 1, will display graphs of data.  


## Example: 

```
mini_analysis(2,11,0.15,29, -8000, 1)
```

This will analyze all waves from sweep2 through sweep11. Each wave will only be analyzed from 0.15 seconds until 29 seconds. The detection cutoff for events is anything with a slope more negative than -8000. When the code is done running it will display analysis graphs.
