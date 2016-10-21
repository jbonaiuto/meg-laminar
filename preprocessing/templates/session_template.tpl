<!-- start session_template.tpl -->
<html>

	<head><title>{PAGETITLE}</title></head>

	<body bgcolor="#ffffff">

		<table border="1" bgcolor="#cccccc" cellpadding="4" cellspacing="0">

			<tr>
    			<td colspan=4><b>{PAGETITLE}</b></td>
			</tr>
            <tr>
                <td colspan=4>Runs</td>
            </tr>
            <tr>
                <td colspan=4>{RUNLIST}</td>
            </tr>
            <tr>
				<td><strong>Instruction - artifact channels</strong></td>
				<td>{INSTRBADCHANNELS}</td>
				<td>Prefiltering<br><img src="{INSTRCHANVARPRESRC}" width="100%"/></td>
				<td>Postfiltering<br><img src="{INSTRCHANVARPOSTSRC}" width="100%"/></td>
			</tr>
			<tr>
				<td><strong>Instruction - artifact trials</strong></td>
				<td>{INSTRBADTRIALS}</td>
				<td>Prefiltering<br><img src="{INSTRTRIALVARPRESRC}" width="100%"/></td>
				<td>Postfiltering<br><img src="{INSTRTRIALVARPOSTSRC}" width="100%"/></td>
			</tr>
			<tr>
				<td><strong>Response - artifact channels</strong></td>
				<td>{RESPBADCHANNELS}</td>
				<td>Prefiltering<br><img src="{RESPCHANVARPRESRC}" width="100%"/></td>
				<td>Postfiltering<br><img src="{RESPCHANVARPOSTSRC}" width="100%"/></td>
			</tr>
			<tr>
				<td><strong>Response - artifact trials</strong></td>
				<td>{RESPBADTRIALS}</td>
				<td>Prefiltering<br><img src="{RESPTRIALVARPRESRC}" width="100%"/></td>
				<td>Postfiltering<br><img src="{RESPTRIALVARPOSTSRC}" width="100%"/></td>
			</tr>
        </table>
	</body>
</html>

<!-- end session_template.tpl -->