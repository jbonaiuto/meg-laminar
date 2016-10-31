<!-- start subject_erp_template.tpl -->
<html>
<head>
<title>{PAGETITLE}</title>
</head>
<body>
<h1>Subject {SUBJECT}</h1>
<h2>{ZEROEVENT}</h2>
WOI: {WOI}
<h2>Mean over sessions (shaded area=std error)</h2>
<table>

	<tr>
		<td><strong>All</strong></td>
		<td><strong>Congruence</strong></td>
		<td><strong>Coherence</strong></td>
		<td><strong>Congruence x Coherence</strong></td>
	</tr>
	<tr>
		<td><img src="{ALLSRC}" width="100%"/></td>
		<td><img src="{CONGRUENCESRC}" width="100%"/></td>
		<td><img src="{COHERENCESRC}" width="100%"/></td>
		<td><img src="{CONGRUENCECOHERENCESRC}" width="100%"/></td>
	</tr>
</table>

<!-- BEGIN session -->
{SESSIONOUT}
<!-- END session -->
<!-- end subject_erp_template.tpl -->
		
</body>
</html>
