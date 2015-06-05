# coding=utf-8
from distancias.items import DistanciasItem
import json
import pkgutil
from scrapy.contrib.spiders import CrawlSpider, Rule
from scrapy.contrib.linkextractors.sgml import SgmlLinkExtractor
from scrapy.selector import Selector

def url_get_ruta(origen, destino):
	return "https://maps.googleapis.com/maps/api/directions/json?origin=" + origen + "&destination=" + destino

class DistanciasSpider(CrawlSpider):
	name = "distancias"
	allowed_domains = ["maps.googleapis.com"]

	url_indices = {}
	start_urls = []
	puntos = pkgutil.get_data('distancias.spiders', 'CH_DU_RM_CL_03.txt').split()
	for i in range(0, len(puntos)):
		for j in range(0, len(puntos)):
			if i != j:
				url = url_get_ruta(puntos[i], puntos[j])
				start_urls.append(url)
				url_indices[url] = [i, j]

	def __init__(self, category=None):
		self.failed_urls = [] # https://stackoverflow.com/questions/13724730/how-to-get-the-scrapy-failure-urls

	def handle_spider_closed(spider, reason):
		with open('failed_urls.txt', 'w') as f:
			f.write(json.dumps(self.failed_urls))

	def parse_start_url(self, response):
		deserializado = json.loads(response.body)

		desde = self.url_indices[response.url][0]
		hasta = self.url_indices[response.url][1]

		if deserializado['status'] == 'OK':
			item = DistanciasItem()
			item['desde'] = desde
			item['hasta'] = hasta
			item['duracion'] = deserializado['routes'][0]['legs'][0]['duration']['value']
			return item
		else:
			self.failed_urls.append([desde, hasta, response.url])
