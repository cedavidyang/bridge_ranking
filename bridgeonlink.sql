BEGIN;
CREATE TABLE network.LA_bridges2 as(
SELECT bridges.*,
	(
	SELECT concat_ws(', ', motorway.gid::int, (motorway.gid+85)::int)
	FROM network.LA_motorway motorway
	ORDER BY ST_Distance(motorway.geom,bridges.geom)
	LIMIT 1) as onlinks
FROM network.LA_bridges bridges);
drop table network.LA_bridges;
alter table network.LA_bridges2 rename to LA_bridges;
COMMIT;
