# -*- coding: utf-8 -*-

# Scrapy settings for distancias project
#
# For simplicity, this file contains only the most important settings by
# default. All the other settings are documented here:
#
#     http://doc.scrapy.org/en/latest/topics/settings.html
#

BOT_NAME = 'distancias'

SPIDER_MODULES = ['distancias.spiders']
NEWSPIDER_MODULE = 'distancias.spiders'

FEED_URI = 'items.csv'
FEED_FORMAT = 'csv'

CSV_DELIMITER = ' '

FIELDS_TO_EXPORT = [
	'desde',
	'hasta',
	'duracion',
]

FEED_EXPORTERS = {
    'csv': 'distancias.exporter.CsvOptionRespectingItemExporter',
}

DOWNLOAD_DELAY = 0.5

# Crawl responsibly by identifying yourself (and your website) on the user-agent
#USER_AGENT = 'distancias (+http://www.yourdomain.com)'
