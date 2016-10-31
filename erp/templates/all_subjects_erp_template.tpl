<!-- start all_subjects_erp_template.tpl -->
<html>
<head>
<title>{PAGETITLE}</title>
</head>
<body>
<h1>{ZEROEVENT}</h1>
WOI: {WOI}
<h2>Mean over all subjects (shaded area=std error)</h2>
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

<!-- BEGIN subject -->
{SUBJECTOUT}
<!-- END subject -->
<!-- end all_subjects_erp_template.tpl -->
</body>
</html>
