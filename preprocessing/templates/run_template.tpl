<!-- start run_template.tpl -->
<html>

	<head><title>{PAGETITLE}</title></head>

	<body bgcolor="#ffffff">

		<table border="1" bgcolor="#cccccc" cellpadding="4" cellspacing="0">

			<tr>
    			<td colspan=4><b>{PAGETITLE}</b></td>
			</tr>
			<tr>
				<td colspan=4><img src="{FIDUCIALSRC}" width="50%"/></td>
			</tr>
			<tr>
				<td colspan=4><img src="{FIDUCIALLASTEVENTSRC}" width="50%"/></td>
			</tr>
			<tr>
				<td colspan=4><img src="{FIDUCIALNOJUMPSRC}" width="50%"/></td>
			</tr>
			<tr>
				<td><strong>Diode channel</strong></td>
				<td>{DIODECHANNEL}</td>
				<td colspan=2></td>
			</tr>
			<tr>
				<td><strong>Diode onsets</strong></td>
				<td>{DIODEONSETS}</td>
				<td colspan=2></td>
			</tr>
            <tr>
				<td><strong>Diode threshold</strong></td>
				<td>{DIODETHRESH}</td>
				<td colspan=2></td>
			</tr>
			<tr>
				<td colspan=4><img src="{DIODESRC}" width="50%"/></td>
			</tr>
			<tr>
				<td><strong>Dots events</strong></td>
				<td>{DOTSEVENTS}</td>
				<td colspan=2></td>
			</tr>
			<tr>
				<td><strong>Instruction events</strong></td>
				<td>{INSTREVENTS}</td>
				<td colspan=2></td>
			</tr>
			<tr>
				<td><strong>Response events</strong></td>
				<td>{RESPEVENTS}</td>
				<td colspan=2></td>
			</tr>
			<tr>
				<td><strong>Channels removed</strong></td>
				<td>{CHANNELSREMOVED}</td>
				<td colspan=2></td>
			</tr>
			<tr>
				<td><strong>Highpass filtered</strong></td>
				<td>{HIGHPASSFREQ}Hz</td>
				<td colspan=2></td>
			</tr>
			<tr>
				<td><strong>Downsampled</strong></td>
				<td>{DOWNSAMPLE}Hz</td>
				<td colspan=2></td>
			</tr>
			<tr>
				<td><strong>Lowpass filtered</strong></td>
				<td>{LOWPASSFREQ}Hz</td>
				<td colspan=2></td>
			</tr>
			<tr>
				<td><strong>Blinks removed</strong></td>
				<td>{BLINKSREMOVED}</td>
				<td colspan=2>
					<img src="{BLINKCMPNTSRC}" width="50%"/>
				</td>
			</tr>
			<tr>
				<td><strong>Instruction-aligned epoch limits</strong></td>
				<td>{INSTREPOCH}</td>
				<td colspan=2></td>
			</tr>
			<tr>
				<td><strong>Response-aligned epoch limits</strong></td>
				<td>{RESPEPOCH}</td>
				<td colspan=2></td>
			</tr>
			<tr>
				<td><strong>No response trials</strong></td>
				<td>{NORESP}</td>
				<td colspan=2></td>
			</tr>
			<tr>
				<td><strong>Incorrect trials</strong></td>
				<td>{INCORRECT}</td>
				<td colspan=2></td>
			</tr>
			

		</table>
	</body>
</html>

<!-- end run_template.tpl -->
