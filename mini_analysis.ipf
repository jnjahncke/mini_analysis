#pragma rtGlobals=1            // Use modern global access method.
#include"justfilters"



Function mini_analysis(startwave,endwave,starttime,endtime, level, plotflag)

	// Input variables:
	variable startwave,endwave,level, plotflag, starttime, endtime
	// Local variables
	Wave W_coef
	Variable/G V_fitOptions=5, V_FitError=0
	variable currentwave=startwave, mini1, mini2, mini3
	variable i, j, srchbefore=0.005, srchafter = 0.01, eventnr=0
	variable tenpercent, ninetypercent, tentime, ninetytime
	Wave decay, iei, Amplitude, RiseTime, InterEventInterval, DecayTime
	string mini_times_wave, mini_diff_wave
	
	Make/O/N=0 Amplitude=NAN
	Make/O/N=0 RiseTime=NAN
	Make/O/N=0 DecayTime=NAN
	Make/O/N=0 InterEventInterval=NAN
	Make/O/N=3 W_coef

	// For rise/decay fit of average wave
	Variable tenpercent_avg, ninetypercent_avg, tentime_avg, ninetytime_avg
	Variable mini_amp_avg, mini_base_avg, basex_avg, peakx_avg, rise_avg, decay_avg, mini_diff_avg
	Variable/G V_fitOptions=5, V_FitError=0
	Wave rise_avgTimeFit, decay_avgTimeFit
	Make /O/N=0 rise_avgTimeFit = NAN, decay_avgTimeFit = NAN

 
 	Make/O/N = (endwave - startwave + 2) minifreq=NAN, miniamp=NAN // Make a new wave with two columns: freq and amp
 	// O = overwright if variable name conflict
 	// N = (rows)
	Make/O/N = (2551, 10000) allevents // Make new wave with 2551 rows, 10000 columns, call it "allevents"
	Make/O/N = 0 allamps = NAN // Make new wave of length 0, call it "allamps"

	for(currentwave = startwave; currentwave <= endwave; currentwave += 1)
		Make/O times
		duplicate/O times base_times // 'times' = sourcewave, 'base_times' = destination wave
		base_times = 0
		duplicate/O $("sweep" + num2str(currentwave)), temp_wave_raw // Make temp_wave_raw to leave UNflitered
		duplicate/O $("sweep" + num2str(currentwave)), temp_wave // Make temp_wave to filter and operate on
		jf_fastfilter(temp_wave, 200, 10) // filter cutoff, poles
		Differentiate temp_wave/D=temp_wave_DIF // first derivative; /D=temp_wave_DIF specifies the name of the wave that will hold the differentiated data
		FindLevels/Q/R=(starttime, endtime)/D=times/EDGE=2 temp_wave_DIF, level // find peaks in first derivative
			// /Q = doesn't abort if no levels are found
			// /R = (starttime, endtime) specifies the X range
			// /D = temp_wave_DIF specifies the name of the destination wave to which to store the level crossing values
			// /EDGE = 2 : searches only for crossing where the Y values are decreasing as the level is crossed
		duplicate/O times, mini_amp, mini_base, basex, peakx, rise, decay, iei, mini_diff, mini_amp_fit, mini_base_fit, basex_fit, peakx_fit, mini_diff_fit // makes new waves from sourcewave 'times'

		
		
		i = 0
		do // cycle through detected minis
			WaveStats/Q/R=(times[i]-srchbefore, times[i]+srchafter) temp_wave_DIF // find the middle of the mini rise
			times[i] = V_minloc
			
			// find mini amplitudes and baseline
			WaveStats/Q/R=(times[i]-srchbefore, times[i]) temp_wave_raw
			mini_base[i] = V_max // y value = baseline value
			basex[i] = V_maxloc // x value = where the baseline occurs in time
			WaveStats/Q/R=(times[i], times[i]+srchafter) temp_wave_raw
			mini_amp[i] = V_min
			peakx[i] = V_minloc
			mini_diff[i] = mini_amp[i] - mini_base[i]
			
			// calculate inter-event-interval
			if (i > 0)
			iei[i] = (basex[i] - basex[i-1])*1000
			else
			iei[i] = NaN
			endif
			
			// find mini amplitudes and baseline FOR KINETICS - idk why but it needs the fited wave
			WaveStats/Q/R=(times[i]-srchbefore, times[i]) temp_wave
			mini_base_fit[i] = V_max // y value = baseline value
			basex_fit[i] = V_maxloc // x value = where the baseline occurs in time
			WaveStats/Q/R=(times[i], times[i]+srchafter) temp_wave
			mini_amp_fit[i] = V_min
			peakx_fit[i] = V_minloc
			mini_diff_fit[i] = mini_amp_fit[i] - mini_base_fit[i]
			
			// find 10% and 90% of amplitude
			tenpercent = 0.1 * mini_diff_fit[i] + mini_base_fit[i]
			ninetypercent = 0.9 * mini_diff_fit[i] + mini_base_fit[i]
			
			// find the elvel crossings for 10% and 90%
			FindLevel/EDGE=0/Q/R=(basex_fit[i], peakx_fit[i]) temp_wave, tenpercent // EDGE=2 looks for decreasing Y values
			tentime = V_LevelX
			FindLevel/EDGE=0/Q/R=(basex_fit[i], peakx_fit[i]) temp_wave, ninetypercent // EDGE=2 looks for decreasing Y values
			ninetytime = V_LevelX
			
			// calculate rise time
			rise[i] = (ninetytime - tentime)*1000
			
			// fit exponential decay
			V_fitOptions=5 // don't update until after the fit is done - same as /N=1
			V_FitError=0 // ignore errors
			CurveFit/Q/N=1/W=2 exp_XOffset temp_wave(peakx_fit[i]+.002, peakx_fit[i]+.02) /D /F={0.95, 4}
			decay[i] = W_coef[2]*1000
			
			// Update Amplitude, RiseTime, and DecayTime lists
			if (decay[i]<50 && decay[i]>0)
			Redimension /N=(numpnts(Amplitude)+1) Amplitude
			Amplitude [numpnts(Amplitude)-1] = mini_diff[i]
			Redimension /N=(numpnts(RiseTime)+1) RiseTime
			RiseTime [numpnts(RiseTime)-1] = rise[i]
			Redimension /N=(numpnts(DecayTime)+1) DecayTime
			DecayTime [numpnts(DecayTime)-1] = decay[i]
			Redimension /N=(numpnts(InterEventInterval)+1) InterEventInterval
			InterEventInterval [numpnts(InterEventInterval)-1] = iei[i]
			// Getting waveform
			duplicate/O/R=(times[i]-0.01, times[i]+0.265) temp_wave_raw, shape
			WaveStats/Q/R=(0,0.265) shape
			shape -= V_avg
			allevents[][eventnr] = shape[p]
			eventnr+=1
			endif
					
			i+=1
			
			while (i < V_LevelsFound)
			
			for (j=1; j<1000; j+=1) // scan the times and look for any duplicate minis
				if (times[j]-times[j-i] <= 0.003)
				DeletePoints j,1,times, mini_amp, mini_base, mini_diff
				endif
				if (decay[i]>50 || decay[i]<0)
				DeletePoints j,1,times, mini_amp, mini_base, mini_diff
				endif
			endfor
			
			for (j=1; j<1000; j+=1) // scan the times again looking for duplicates
				if(times[j]-times[j-i] <=0.003)
				DeletePoints j,1,times, mini_amp, mini_base, mini_diff
				endif
			endfor


		WaveStats/Q mini_diff // calculate average mini frequency for each sweep
		minifreq[currentwave-startwave] = V_npnts/(endtime-starttime)
		miniamp[currentwave-startwave]=V_avg
		
		if(plotflag == 1 && currentwave==endwave) // show the resulting waves and points
			display /K=1/W=(5,20,315,194) temp_wave_raw
			SetAxis bottom starttime,endtime
			appendtograph mini_amp vs times
			appendtograph mini_base vs times
			ModifyGraph mode(mini_amp)=3, marker(mini_amp)=19
			ModifyGraph useMrkStrokeRGB(mini_amp)=1
			ModifyGraph rgb(mini_amp) = (0,0,0)
			ModifyGraph mode(mini_base)=3, marker(mini_base)=9, msize(mini_base)=3
		
			
			Edit/K=1 'minifreq';DelayUpdate // Creates the frequency table
			AppendToTable 'miniamp';DelayUpdate // adds amplitude column to the top table
		endif

      
      // Make new waves for mini_times (tsweep) and mini_diff (asweep)
      mini_times_wave="t"+num2str(currentwave)
      mini_diff_wave="a"+num2str(currentwave)
      make/O/N=0 $mini_times_wave
      make/O/N=0 $mini_diff_wave
      wave app1 = $mini_times_wave
      wave app2 = $mini_diff_wave
      duplicate/O times, app1
      duplicate/O mini_diff, app2                                   
      
      concatenate {app2}, allamps  
       
	endfor // end cycling through waves

	WaveStats/Q minifreq // calculate average mini frequency for all sweeps
	minifreq[endwave-startwave+1]=V_avg
	WaveStats/Q miniamp // calculate average mini amplitude for all sweeps
	miniamp[endwave-startwave+1]=V_avg              

	duplicate/O/R=[][0,eventnr] allevents, events
	killwaves allevents
	SetScale/P x 0,10e-05, "s", events        

	// showing all events
	i=1
	display/K=1 events
	duplicate/O/R=[][0] events,avg_event
	do
		appendtograph events[][i]
		avg_event = (avg_event*i+events[p][i])/(i+1)
		i+=1
	while(i<=eventnr)
	redimension/N=-1 avg_event
	appendtograph avg_event
	ModifyGraph rgb(avg_event)=(0,0,0)
	SetAxis bottom *,0.075
	
	// scatter plot amplitude vs rise time
	Display/K=1 RiseTime vs Amplitude
	ModifyGraph mode=2,rgb=(0,0,0)
	SetAxis bottom 0,-800
	ModifyGraph lsize=2

	if(plotflag == 1)
			Make/o/N=1 histResult  //make a histogram of all event amplitudes
		Allamps *=-1
		Histogram/B={-100,10,310} allamps, histResult // edited 3/23/2020 to go out to 3000pA
		Display/K=1 histResult; ModifyGraph mode=5
		ModifyGraph log(left)=1
	endif

	// show mini amplitude vs rise time
	Edit /k=1 Amplitude,RiseTime,DecayTime,InterEventInterval
	
	// make a histogram of inter event intervals
	if(plotflag == 1)
			Make/o/N=1 histiei
		Histogram /B={0.000001,2,500} InterEventInterval, histiei
		Display/K=1 histiei; ModifyGraph mode=5
	endif
	
	// calculate rise and decay kinetics of the average event
	// find base, peak
	WaveStats/Q/R=(0.01-srchbefore, 0.01) avg_event
	mini_base_avg = V_max
	basex_avg = V_maxloc
	
	WaveStats/Q/R=(0.01, 0.01+srchafter) avg_event
	mini_amp_avg = V_min
	peakx_avg = V_minloc
	mini_diff_avg = mini_amp_avg - mini_base_avg
	
	// find 10% and 90%
	tenpercent_avg = 0.1 * mini_diff_avg + mini_base_avg
	ninetypercent_avg = 0.9 * mini_diff_avg + mini_base_avg
	
	FindLevel/EDGE=0/Q/R=(basex_avg, peakx_avg) avg_event, tenpercent_avg
	tentime_avg = V_LevelX
	FindLevel/EDGE=0/Q/R=(basex_avg, peakx_avg) avg_event, ninetypercent_avg
	ninetytime_avg = V_LevelX
	
	// calculate rise_avg time
	rise_avg = (ninetytime_avg - tentime_avg)*1000
	
	// calculate decay_avg time
	V_fitOptions=5
	V_FitError=0
	CurveFit/Q/N=1/W=2 exp_Xoffset avg_event(peakx_avg+0.002, peakx_avg+0.02) /D /F={0.95, 4}
	decay_avg = W_coef[2]*1000
	
	// update table
	Redimension /N=(numpnts(rise_avgTimeFit)+1) rise_avgTimeFit
	rise_avgTimeFit [numpnts(rise_avgTimeFit)-1] = rise_avg
	Redimension /N=(numpnts(decay_avgTimeFit)+1) decay_avgTimeFit
	decay_avgTimeFit [numpnts(decay_avgTimeFit)-1] = decay_avg
	
	// show rise_avg and decay_avg
	Edit /k=1 rise_avgTimeFit,decay_avgTimeFit
 	
End
