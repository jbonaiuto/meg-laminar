<!-- start session_erp_template.tpl -->

<h2>Session {SESSIONNUM}</h2>
<h3>Mean over trials (shaded area=std error)</h3>
<table>
    <tr>

		<td><strong>Channel</strong>: {CHAN}</td>
		<td><img src="{CHANMAP}" width="50%"/></td>
		<td colspan=2></td>

	</tr>
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

<!-- end session_erp_template.tpl -->
		
