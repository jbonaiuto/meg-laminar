<!-- start all_subjects_sensor_tf_template.tpl -->
<html>

	<head><title>{PAGETITLE}</title></head>

	<body bgcolor="#ffffff">
		<h1>{PAGETITLE}</h1>
		<h1>Events</h1>
		<ul>
			<li><a href="#dots_aligned">Dots-aligned</a></li>
			<li><a href="#instruction_aligned">Instruction-aligned</a></li>
			<li><a href="#response_aligned">Response-aligned</a></li>
		</ul>
		<h1>Subjects</h1>
		{SUBJLIST}
		<a name="dots_aligned"><h2>Dots-aligned</h2></a>
		<table border="1" bgcolor="#cccccc" cellpadding="4" cellspacing="0">

			<tr>
				<td><b>Frequency Range</b></td>
				<td><b>F-test</b></td>
    				<td><b>Positive</b></td>
				<td><b>Negative</b></td>
			</tr>
			<tr>
				<td></td>
				<td><img src="{DOTSFTFSRC}" width="100%"/></td>
				<td><img src="{DOTSPOSITIVETFSRC}" width="100%"/></td>
				<td><img src="{DOTSNEGATIVETFSRC}" width="100%"/></td>
			</tr>	
			<tr>
				<td>Broadband</td>
				<td><img src="{DOTSFSFBROADSRC}" width="100%"/></td>
				<td><img src="{DOTSPOSITIVESFBROADSRC}" width="100%"/></td>
				<td><img src="{DOTSNEGATIVESFBROADSRC}" width="100%"/></td>
            		<tr>
			<tr>
				<td>Alpha</td>
				<td><img src="{DOTSFSFALPHASRC}" width="100%"/></td>
				<td><img src="{DOTSPOSITIVESFALPHASRC}" width="100%"/></td>
				<td><img src="{DOTSNEGATIVESFALPHASRC}" width="100%"/></td>
            		<tr>
			<tr>
				<td>Beta</td>
				<td><img src="{DOTSFSFBETASRC}" width="100%"/></td>
				<td><img src="{DOTSPOSITIVESFBETASRC}" width="100%"/></td>
				<td><img src="{DOTSNEGATIVESFBETASRC}" width="100%"/></td>
            		<tr>
			<tr>
				<td>Gamma</td>
				<td><img src="{DOTSFSFGAMMASRC}" width="100%"/></td>
				<td><img src="{DOTSPOSITIVESFGAMMASRC}" width="100%"/></td>
				<td><img src="{DOTSNEGATIVESFGAMMASRC}" width="100%"/></td>
            		<tr>
		</table>
		<a name="instruction_aligned"><h2>Instruction-aligned</h2></a>
		<table border="1" bgcolor="#cccccc" cellpadding="4" cellspacing="0">

			<tr>
				<td><b>Frequency Range</b></td>
				<td><b>F-test</b></td>
    				<td><b>Positive</b></td>
				<td><b>Negative</b></td>
			</tr>
			<tr>
				<td></td>
				<td><img src="{INSTRFTFSRC}" width="100%"/></td>
				<td><img src="{INSTRPOSITIVETFSRC}" width="100%"/></td>
				<td><img src="{INSTRNEGATIVETFSRC}" width="100%"/></td>
			</tr>	
			<tr>
				<td>Broadband</td>
				<td><img src="{INSTRFSFBROADSRC}" width="100%"/></td>
				<td><img src="{INSTRPOSITIVESFBROADSRC}" width="100%"/></td>
				<td><img src="{INSTRNEGATIVESFBROADSRC}" width="100%"/></td>
            		<tr>
			<tr>
				<td>Alpha</td>
				<td><img src="{INSTRFSFALPHASRC}" width="100%"/></td>
				<td><img src="{INSTRPOSITIVESFALPHASRC}" width="100%"/></td>
				<td><img src="{INSTRNEGATIVESFALPHASRC}" width="100%"/></td>
            		<tr>
			<tr>
				<td>Beta</td>
				<td><img src="{INSTRFSFBETASRC}" width="100%"/></td>
				<td><img src="{INSTRPOSITIVESFBETASRC}" width="100%"/></td>
				<td><img src="{INSTRNEGATIVESFBETASRC}" width="100%"/></td>
            		<tr>
			<tr>
				<td>Gamma</td>
				<td><img src="{INSTRFSFGAMMASRC}" width="100%"/></td>
				<td><img src="{INSTRPOSITIVESFGAMMASRC}" width="100%"/></td>
				<td><img src="{INSTRNEGATIVESFGAMMASRC}" width="100%"/></td>
            		<tr>
		</table>
		<a name="response_aligned"><h2>Response-aligned</h2></a>
		<table border="1" bgcolor="#cccccc" cellpadding="4" cellspacing="0">

			<tr>
				<td><b>Frequency Range</b></td>
				<td><b>F-test</b></td>
    				<td><b>Positive</b></td>
				<td><b>Negative</b></td>
			</tr>
			<tr>
				<td></td>
				<td><img src="{RESPFTFSRC}" width="100%"/></td>
				<td><img src="{RESPPOSITIVETFSRC}" width="100%"/></td>
				<td><img src="{RESPNEGATIVETFSRC}" width="100%"/></td>
			</tr>	
			<tr>
				<td>Broadband</td>
				<td><img src="{RESPFSFBROADSRC}" width="100%"/></td>
				<td><img src="{RESPPOSITIVESFBROADSRC}" width="100%"/></td>
				<td><img src="{RESPNEGATIVESFBROADSRC}" width="100%"/></td>
            		<tr>
			<tr>
				<td>Alpha</td>
				<td><img src="{RESPFSFALPHASRC}" width="100%"/></td>
				<td><img src="{RESPPOSITIVESFALPHASRC}" width="100%"/></td>
				<td><img src="{RESPNEGATIVESFALPHASRC}" width="100%"/></td>
            		<tr>
			<tr>
				<td>Beta</td>
				<td><img src="{RESPFSFBETASRC}" width="100%"/></td>
				<td><img src="{RESPPOSITIVESFBETASRC}" width="100%"/></td>
				<td><img src="{RESPNEGATIVESFBETASRC}" width="100%"/></td>
            		<tr>
			<tr>
				<td>Gamma</td>
				<td><img src="{RESPFSFGAMMASRC}" width="100%"/></td>
				<td><img src="{RESPPOSITIVESFGAMMASRC}" width="100%"/></td>
				<td><img src="{RESPNEGATIVESFGAMMASRC}" width="100%"/></td>
            		<tr>
		</table>
	</body>
</html>

<!-- end all_subjects_sensor_tf_template.tpl -->
