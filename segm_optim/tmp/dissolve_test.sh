ogr2ogr output2.shp park_segm_2017_1.shp -dialect sqlite -sql "SELECT ST_Union(geometry), segm_id FROM park_segm_2017_1 GROUP BY segm_id"


saga_cmd shapes_polygons 5 -POLYGONS 'park_segm_2017_1.shp' -FIELD_1 'segm_id' -DISSOLVED 'output.shp'


