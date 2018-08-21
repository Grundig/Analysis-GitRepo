url = 'https://www.yr.no/place/Norway/Hordaland/Bergen/Bergen/hour_by_hour_detailed.html';
html = urlread(url);
a = strfind(html,'<td title="Pressure:');
pressure = html(a(1)+21:a(1)+26);